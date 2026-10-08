"""Equal-total-nt experiments. See docs/BENCHMARK.md for statistical scope."""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import re
import subprocess
import sys
import time
import zipfile
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import edlib
import numpy as np
import yaml

from .channel import ChannelParams, simulate, copy_weights
from .codecs import make_codec
from .config import deep_merge, load_config
from .container import Header, prepare
from .decoder import decode_reads
from .ecc import CRC_BITS, plan_layout
from .encoder import encode_bytes, write_fasta
from .metrics import pool_metrics
from .oligo import HEADER_INDEX_SPACE, OligoGeometry, header_fragments
from .primers import load_primers
from .screening import ScreeningConfig, screen_pool

ROOT = Path(__file__).resolve().parent.parent


def wilson(successes: int, n: int) -> tuple[float, float]:
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("Wilson interval requires 0 <= successes <= n and n > 0")
    z = 1.959963984540054
    p, den = successes / n, 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0., mid - half), min(1., mid + half)


def derived_seed(root: int, *keys) -> int:
    # Stable labels: adding/reordering arms never changes the other arms' RNG streams.
    digest = hashlib.sha256(json.dumps(keys, separators=(",", ":")).encode()).digest()
    words = np.frombuffer(digest, dtype="<u4").tolist()
    return int(np.random.SeedSequence([root, *words]).generate_state(1)[0])


def _json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + "\n")
    tmp.replace(path)


def _csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    tmp = path.with_suffix(".tmp")
    with tmp.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _shape(data, primers, cfg):
    """Plan header size and capacity without running an expensive steering search."""
    stored, info = prepare(data, **cfg["container"])
    geom = OligoGeometry(cfg["oligo_len"], len(primers.forward), len(primers.reverse), cfg.get("seed_trits", 3))
    geom.validate()
    codec = make_codec(cfg["codec"]["name"], cfg["codec"].get("params"))
    codec.bind(primers, cfg["screening"])
    bits = codec.frame_bits(geom.payload_nt) - CRC_BITS
    if bits < 8:
        raise ValueError("codec capacity must be at least 8 payload bits")
    if getattr(codec, "rateless", False):
        codec.seg_bytes = codec.seg_bytes_for(geom.payload_nt)
        codec.k = math.ceil(len(stored) / codec.seg_bytes)
        bits = codec.seg_bytes * 8
    h = Header(info.file_type, info.original_len, info.stored_len, info.compressed,
               info.crc32, cfg["global_seed"], codec.codec_id, bits, geom.body_nt, geom.seed_trits,
               1, 0, codec.params())
    return stored, geom, codec, h


def _header_count(h, geom, copies):
    n = len(header_fragments(h.pack(), geom.header_payload_nt)) * copies
    if not 1 <= n <= HEADER_INDEX_SPACE:
        raise ValueError("header copies exceed the reserved index space")
    return n


def match_budget(data, primers, base, arms, max_parity_permille=2000):
    """Smallest exactly shared pool size reachable by all RS arms, then LT matching.

    Search integer parity ratios, retaining the lowest ratio for each pool size.
    Header fragments, primers, bit packing and row padding all count. No filler oligos.
    """
    choices, fountains = {}, {}
    for arm in arms:
        cfg = deep_merge(base, {"codec": arm["codec"]})
        stored, geom, codec, h = _shape(data, primers, cfg)
        if not stored:
            raise ValueError("benchmarks require nonempty stored payloads")
        if getattr(codec, "rateless", False):
            fountains[arm["name"]] = (cfg, geom, codec, h)
            continue
        counts = {}
        for perm in range(base["ecc"]["parity_permille"], max_parity_permille + 1):
            try:
                layout = plan_layout(len(stored), h.payload_bits, perm, base["ecc"].get("k_max"))
            except ValueError:
                continue
            header = replace(h, rs_k_max=layout.k_max, rs_parity_permille=perm)
            n = layout.n_oligos + _header_count(header, geom, cfg["header_copies"])
            counts.setdefault(n, deep_merge(cfg, {"ecc": {"parity_permille": perm, "k_max": layout.k_max}}))
        choices[arm["name"]] = counts
    if not choices:
        raise ValueError("at least one RS arm is required for budget matching")
    shared = set.intersection(*(set(x) for x in choices.values()))
    for target in sorted(shared):
        selected = {name: opts[target] for name, opts in choices.items()}
        for name, (cfg, geom, codec, h) in fountains.items():
            # Floating-point overhead changes JSON length, hence header count. Solve both.
            found = None
            for headers in range(1, HEADER_INDEX_SPACE + 1):
                droplets = target - headers
                if droplets < math.ceil(codec.k * (1 + codec.overhead)):
                    continue
                overhead = (droplets - 0.5) / codec.k - 1
                params = {**cfg["codec"].get("params", {}), "overhead": overhead}
                candidate_h = replace(h, codec_params={**h.codec_params, "overhead": overhead})
                if _header_count(candidate_h, geom, cfg["header_copies"]) == headers:
                    found = deep_merge(cfg, {"codec": {"params": params}})
                    break
            if found is None:
                break
            selected[name] = found
        else:
            return target, selected
    raise ValueError("no exact common nucleotide budget; increase max_parity_permille or adjust arms/payload size")


def oligo_metrics(result, enc, primers):
    errors, present = [], []
    n_data = len(enc.oligos) - enc.n_header_oligos
    for idx, seq in enc.oligos:
        if idx < HEADER_INDEX_SPACE:
            continue
        body = seq[len(primers.forward):-len(primers.tail)]
        got = result.consensus.get(idx)
        if got is None:
            errors.append(1.)
        else:
            error = edlib.align(got, body, task="distance")["editDistance"] / len(body)
            errors.append(error)
            present.append(error)
    valid = result.report.get("data_oligos_valid", result.report.get("droplets_valid", 0))
    return {"body_error_mean": float(np.mean(errors)) if errors else 0.,
            "body_error_present_mean": float(np.mean(present)) if present else None,
            "missing_oligos": n_data - len(present), "valid_oligos": valid,
            "unrecovered_oligo_fraction": 1 - valid / n_data if n_data else 0.}


def summarize(rows):
    groups = defaultdict(list)
    keys = ("payload", "payload_seed", "arm", "channel", "strength", "depth")
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    out = []
    for key, values in sorted(groups.items()):
        n, wins = len(values), sum(v["ok"] for v in values)
        lo, hi = wilson(wins, n)
        out.append({**dict(zip(keys, key)), "n": n, "recovered": wins, "recovery_fraction": wins / n,
                    "wilson_low": lo, "wilson_high": hi,
                    **{k: float(np.mean([v[k] for v in values])) for k in
                       ("body_error_mean", "missing_oligos", "unrecovered_oligo_fraction", "actual_depth")}})
    return out


def depth_boundary(probe, max_depth):
    """Bisection for observed all-seed recovery. Returns None when right-censored.

    Recovery is stochastic; this is an empirical boundary, not a guaranteed minimum.
    The caller saves every probe and its interval, including nonmonotone observations.
    """
    if not probe(max_depth):
        return None
    lo, hi = 0, max_depth
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if probe(mid):
            hi = mid
        else:
            lo = mid
    return hi


def _reads(enc, channel, depth, seed, work):
    seqs = [s for _, s in enc.oligos]
    params = dict(channel.get("params", {}))
    if channel["kind"] == "ids":
        cp = ChannelParams.from_dict({**params, "mean_coverage": depth, "seed": seed})
        return [r.seq for r in simulate(seqs, cp, np.random.default_rng(seed))[0]]
    from . import nanopore
    identity = params.pop("identity", nanopore.BADREAD_IDENTITY)
    cp = ChannelParams.from_dict({**params, "mean_coverage": depth})
    rng = np.random.default_rng(seed)
    w = copy_weights(seqs, cp, rng)
    w[rng.random(len(w)) < cp.dropout] = 0
    if w.sum() == 0:
        return []
    copies = w * (depth * len(w) / w.sum())
    reads, _ = nanopore.badread_reads(seqs, copies, work, seed, identity=identity)
    return [seq for _, seq in reads]


def _validate(spec):
    allowed = {"mode", "seed", "encode_config", "encode", "payloads", "payload_seeds", "channel_seeds",
               "arms", "channels", "depths", "max_depth", "minimum_depth", "max_parity_permille", "plots"}
    if not isinstance(spec, dict) or set(spec) - allowed:
        raise ValueError("benchmark must be a mapping with known keys (see docs/BENCHMARK.md)")
    if spec.get("mode", "study") not in ("study", "smoke"):
        raise ValueError("mode must be study or smoke")
    n = spec.get("channel_seeds", 20)
    if not isinstance(n, int) or n < (1 if spec.get("mode") == "smoke" else 20):
        raise ValueError("study mode requires >=20 channel_seeds; smoke requires >=1")
    for key in ("payloads", "arms", "channels"):
        items = spec.get(key)
        if not isinstance(items, list) or not items:
            raise ValueError(f"{key} must be a nonempty list")
        names = [x.get("name", "") for x in items]
        if len(set(names)) != len(names) or any(not re.fullmatch(r"[A-Za-z0-9_-]+", n) for n in names):
            raise ValueError(f"{key} names must be unique letters/digits/hyphens/underscores")
    for p in spec["payloads"]:
        if ("path" in p) == ("bytes" in p) or ("bytes" in p and (not isinstance(p["bytes"], int) or p["bytes"] < 1)):
            raise ValueError("each payload needs either a path or positive bytes")
    seeds = spec.get("payload_seeds", [0])
    if not seeds or len(set(seeds)) != len(seeds) or any(not isinstance(s, int) or s < 0 for s in seeds):
        raise ValueError("payload_seeds must be unique nonnegative integers")
    if not isinstance(spec.get("seed", 0), int) or spec.get("seed", 0) < 0:
        raise ValueError("seed must be a nonnegative integer")
    if not isinstance(spec.get("max_parity_permille", 2000), int) or not 0 <= spec.get("max_parity_permille", 2000) <= 65535:
        raise ValueError("max_parity_permille must be an integer in [0, 65535]")
    for arm in spec["arms"]:
        if not isinstance(arm.get("codec"), dict) or "name" not in arm["codec"]:
            raise ValueError("each arm must provide a codec mapping with a name")
        make_codec(arm["codec"]["name"], arm["codec"].get("params"))
    depths = spec.get("depths", [3, 5, 10])
    maximum = spec.get("max_depth", max(depths) if depths else 0)
    if not depths or any(not isinstance(d, int) or d < 1 for d in depths) or len(set(depths)) != len(depths):
        raise ValueError("depths must be unique positive integers")
    if not isinstance(maximum, int) or maximum < max(depths):
        raise ValueError("max_depth must be an integer >= all requested depths")
    for c in spec["channels"]:
        if c["kind"] not in ("ids", "badread"):
            raise ValueError("benchmark channels must be ids or badread; use nanopore_e2e.py for signal stress tests")
        params = dict(c.get("params", {}))
        params.pop("identity", None)
        cp = ChannelParams.from_dict(params)
        if cp.conf_weighted or cp.coverage == "model":
            raise ValueError("main benchmark excludes circular and assumed accessibility channels; use sensitivity scripts")
        if c["kind"] == "badread" and (cp.p_sub or cp.p_ins or cp.p_del or cp.rc_frac != 0.5):
            raise ValueError("Badread sets its own errors and 50% strand mixture; IDS rates/rc_frac are unsupported")
    return spec


def _plots(rows, out, mode):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    groups = defaultdict(list)
    for r in summarize(rows):
        groups[(r["payload"], r["payload_seed"], r["channel"], r["strength"])].append(r)
    for key, group in groups.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        arms = sorted({r["arm"] for r in group})
        single_depth = len({r["depth"] for r in group}) == 1
        for position, arm in enumerate(arms):
            values = sorted((r for r in group if r["arm"] == arm), key=lambda r: r["depth"])
            y = np.array([r["recovery_fraction"] for r in values])
            bounds = np.array([[r["wilson_low"], r["wilson_high"]] for r in values])
            x = [position] if single_depth else [r["depth"] for r in values]
            ax.errorbar(x, y, yerr=np.maximum(0, np.array([y-bounds[:, 0], bounds[:, 1]-y])),
                        marker="o", capsize=3, label=arm)
        ax.set(xlabel=f"Codec (requested depth {group[0]['depth']})" if single_depth else "Requested mean sequencing depth",
               ylabel="File recovery fraction (Wilson 95% CI)",
               ylim=(-.05, 1.05), title=f"{mode.upper()} | {' / '.join(map(str, key))}")
        if single_depth:
            ax.set_xticks(range(len(arms)), arms, rotation=20, ha="right")
        else:
            ax.legend(fontsize=8)
        ax.grid(alpha=.2)
        fig.tight_layout()
        stem = "recovery_" + "_".join(map(str, key))
        for ext in ("png", "svg"):
            fig.savefig(out / f"{stem}.{ext}", dpi=160)
        plt.close(fig)


def run_benchmark(config_path, out_dir, progress=print):
    path, out = Path(config_path).resolve(), Path(out_dir).resolve()
    spec = _validate(yaml.safe_load(path.read_text()))
    base_path = path.parent / spec["encode_config"] if spec.get("encode_config") else None
    base = load_config(base_path, spec.get("encode", {}))
    primers = load_primers(base["primers"])
    # Fail before spending time encoding if an optional channel tool is absent.
    if any(c["kind"] == "badread" for c in spec["channels"]):
        from . import nanopore
        if not nanopore.BADREAD.is_file():
            raise ValueError("Badread is not installed; see docs/BENCHMARK.md")
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"output directory is not empty: {out}; choose a new directory")
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.yaml").write_text(path.read_text())
    _json(out / "encode_config.json", base)
    root_seed = spec.get("seed", 0)
    provenance = {"status": "running", "mode": spec.get("mode", "study"), "python": sys.version,
                  "seed": root_seed, "config_sha256": _sha(path), "primer_sha256": _sha(Path(base["primers"])),
                  "packages": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()},
                  "source_sha256": {str(p.relative_to(ROOT)): _sha(p) for p in sorted((ROOT / "dnastore").rglob("*.py"))},
                  "payloads": [], "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    provenance["git_sha"] = git.stdout.strip() if git.returncode == 0 else None
    model_prov = ROOT / "data/kmer_models/PROVENANCE.json"
    if model_prov.exists():
        provenance["kmer_models"] = json.loads(model_prov.read_text())
    # Source hashes detect edits; the archive also preserves the code when Git is absent.
    with zipfile.ZipFile(out / "source.zip", "w", compression=zipfile.ZIP_DEFLATED) as archive:
        sources = list((ROOT / "dnastore").rglob("*.py")) + list((ROOT / "dnastore/resources").glob("*.yaml"))
        sources += [p for p in (ROOT / "pyproject.toml", ROOT / "README.md") if p.exists()]
        for source in sorted(sources):
            archive.write(source, str(source.relative_to(ROOT)))
    provenance["source_archive_sha256"] = _sha(out / "source.zip")
    _json(out / "provenance.json", provenance)
    rows, encoded, boundaries = [], [], []
    try:
        for payload in spec["payloads"]:
            for payload_seed in spec.get("payload_seeds", [0]):
                tag = f"{payload['name']}_{payload_seed}"
                data_seed = derived_seed(root_seed, "payload", payload["name"], payload_seed)
                data = ((path.parent / payload["path"]).read_bytes() if "path" in payload
                        else np.random.default_rng(data_seed).bytes(payload["bytes"]))
                provenance["payloads"].append({"name": tag, "bytes": len(data), "seed": data_seed,
                                                "sha256": hashlib.sha256(data).hexdigest(), "source": payload})
                cfg = deep_merge(base, {"global_seed": derived_seed(root_seed, "encoder", tag)})
                target, configs = match_budget(data, primers, cfg, spec["arms"], spec.get("max_parity_permille", 2000))
                progress(f"{tag}: matched {target} oligos / {target * cfg['oligo_len']} nt per arm")
                for arm in spec["arms"]:
                    name = arm["name"]
                    arm_dir = out / "pools" / tag / name
                    arm_dir.mkdir(parents=True)
                    enc = encode_bytes(data, primers, configs[name])
                    if len(enc.oligos) != target or any(len(s) != cfg["oligo_len"] for _, s in enc.oligos):
                        raise RuntimeError("encoded pool does not match planned nucleotide budget")
                    write_fasta(enc.oligos, str(arm_dir / "oligos.fasta"))
                    clean = decode_reads([s for _, s in enc.oligos], primers, cfg["oligo_len"],
                                         workers=cfg["workers"], strength="erasure")
                    screens = screen_pool([s for _, s in enc.oligos], primers, ScreeningConfig.from_dict(cfg["screening"]),
                                          cfg["compute_mfe"], cfg["workers"])
                    metrics = pool_metrics(enc)
                    encoded.append({"payload": payload["name"], "payload_seed": payload_seed, "arm": name,
                                    "noiseless_ok": clean.ok and clean.data == data,
                                    "violating_oligos": sum(bool(s["violations"]) for s in screens), **metrics})
                    _json(arm_dir / "manifest.json", {"config": configs[name], "metrics": metrics,
                            "noiseless_report": clean.report, "screening": screens, "encoder_info": enc.payload_info,
                            "fasta_sha256": _sha(arm_dir / "oligos.fasta")})
                    _csv(out / "encode.csv", encoded)
                    progress(f"  {name}: noiseless recovery={clean.ok and clean.data == data}")
                    for channel in spec["channels"]:
                        cache = {}
                        def probe(depth):
                            if depth in cache:
                                return cache[depth]
                            success = {"erasure": True, "repair": True}
                            for rep in range(spec.get("channel_seeds", 20)):
                                seed = derived_seed(root_seed, "channel", tag, channel["name"], depth, rep)
                                work = out / "channels" / tag / name / channel["name"] / f"d{depth}_s{rep}"
                                reads = _reads(enc, channel, depth, seed, work)
                                for strength in ("erasure", "repair"):
                                    start = time.perf_counter()
                                    result = decode_reads(reads, primers, cfg["oligo_len"], workers=cfg["workers"], strength=strength)
                                    ok = result.ok and result.data == data
                                    success[strength] &= ok
                                    row = {"payload": payload["name"], "payload_seed": payload_seed, "arm": name,
                                           "channel": channel["name"], "strength": strength, "depth": depth, "rep": rep,
                                           "seed": seed, "ok": ok, "reads": len(reads), "actual_depth": len(reads) / target,
                                           "total_nt": metrics["total_nt"], "decode_s": time.perf_counter() - start,
                                           "error": result.report.get("error", ""), **oligo_metrics(result, enc, primers)}
                                    rows.append(row)
                            cache[depth] = success
                            _csv(out / "trials.csv", rows)
                            _csv(out / "summary.csv", summarize(rows))
                            progress(f"    {channel['name']} depth={depth}: " + ", ".join(
                                f"{st}={sum(r['ok'] for r in rows[-2 * spec.get('channel_seeds', 20):] if r['strength'] == st)}/{spec.get('channel_seeds', 20)}"
                                for st in success))
                            return success
                        for depth in sorted(spec.get("depths", [3, 5, 10])):
                            probe(depth)
                        if spec.get("minimum_depth", True):
                            maximum = spec.get("max_depth", max(spec.get("depths", [3, 5, 10])))
                            found = {st: depth_boundary(lambda d, st=st: probe(d)[st], maximum)
                                     for st in ("erasure", "repair")}
                            for st in ("erasure", "repair"):
                                boundary = found[st]
                                observed = [d for d, v in cache.items() if v[st]]
                                nonmonotone = any(cache[a][st] and not cache[b][st] for a in cache for b in cache if a < b)
                                boundaries.append({"payload": payload["name"], "payload_seed": payload_seed, "arm": name,
                                    "channel": channel["name"], "strength": st, "bisection_depth": boundary,
                                    "lowest_observed_all_success_depth": min(observed) if observed else None,
                                    "max_depth": maximum, "right_censored": boundary is None, "nonmonotone": nonmonotone,
                                    "channel_seeds": spec.get("channel_seeds", 20)})
                            _csv(out / "minimum_depth.csv", boundaries)
        if spec.get("plots", True):
            _plots(rows, out, provenance["mode"])
        provenance["status"] = "complete"
        provenance["trials"] = len(rows)
    except BaseException as exc:
        provenance["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        provenance["error"] = str(exc)
        raise
    finally:
        provenance["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        _json(out / "provenance.json", provenance)
    return rows
