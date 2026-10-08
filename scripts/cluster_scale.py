"""Clustering scale check at Phase 4 file sizes (thousands of oligos per pool).

Per arm (one ``--bytes`` random payload, seed 0):
1. Encode and save the pool (``pools/<arm>.fa``, reused by later phases).
2. **Distinct-oligo distance:** exact minimum edit distance between oligo bodies over ALL
   pairs (edlib with a bound of ``--bound``·L, so far pairs are cheap), plus the number of
   pairs below 0.30 / 0.33 / 0.36 · L. Phase 3's 0.36 came from ~245 oligos.
3. **Clustering on reads** for each channel (``badread``: calibrated R10.4.1 read model;
   ``ids9``: uniform IDS at 9% total, i.e. the squigulator raw error, as a stress test), at
   lognormal (sigma 0.27) mean depth ``--depth``, with the decoder's own orientation and
   :func:`kmer_cluster` call. Measured against the truth:
   - merge rate: fraction of clusters (>= 2 reads) holding reads of more than one oligo, and
     the fraction of reads that are a minority oligo in their cluster;
   - split rate: fraction of oligos whose reads fall in >= 2 clusters, and clusters per oligo.
4. **Full decode** (``repair`` and ``erasure`` strengths): file recovery, and
   ``groups_multi_index``, the number of clusters that yielded more than one CRC-valid oligo
   (a merge that the per-read-index sub-grouping split back apart).

Usage: .venv/bin/python scripts/cluster_scale.py --bytes 100000 --out results/phase3/cluster_scale
"""

from __future__ import annotations

import argparse
import json
import random
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import edlib
import numpy as np
import pandas as pd

from dnastore import nanopore
from dnastore.channel import ChannelParams, simulate
from dnastore.config import load_config
from dnastore.consensus import kmer_cluster, locate_and_orient
from dnastore.decoder import decode_reads
from dnastore.encoder import encode_bytes, write_fasta
from dnastore.primers import load_primers

ROOT = Path(__file__).resolve().parent.parent
ARMS = [
    ("naive2bit-whiten", "naive2bit", {"whiten": True}),
    ("goldman", "goldman", {}),
    ("steering-A-P8", "steering", {"P": 8, "alphabet": "A"}),
    ("steering-C-P8", "steering", {"P": 8, "alphabet": "C"}),
    ("steering-B-P6", "steering", {"P": 6, "alphabet": "B"}),
    ("fountain", "fountain", {"overhead": 0.15}),
]
_B: list[str] = []


def _init(bodies: list[str]) -> None:
    _B[:] = bodies


def _row_min(args: tuple[int, int]) -> tuple[int, list[int]]:
    i, bound = args
    a, best, below = _B[i], 10**9, []
    for j in range(i + 1, len(_B)):
        r = edlib.align(a, _B[j], task="distance", k=bound)["editDistance"]
        if r >= 0:
            below.append(r)
            best = min(best, r)
    return best, below


def min_distance(bodies: list[str], bound_frac: float, workers: int) -> dict:
    L = len(bodies[0])
    bound = int(bound_frac * L)
    with ProcessPoolExecutor(workers, initializer=_init, initargs=(bodies,)) as ex:
        res = list(ex.map(_row_min, [(i, bound) for i in range(len(bodies))], chunksize=16))
    allb = [d for _, b in res for d in b]
    m = min((x for x, _ in res), default=10**9)
    return {"n_bodies": len(bodies), "pairs": len(bodies) * (len(bodies) - 1) // 2,
            "min_ed": m if m < 10**9 else f">{bound}", "min_ed_frac": m / L if m < 10**9 else None,
            **{f"pairs_le_{t}": sum(d <= t * L for d in allb) for t in (0.30, 0.33, 0.36)}}


def channel_reads(label: str, oligos: list[str], depth: float, seed: int, out: Path):
    rng = np.random.default_rng(seed)
    copies = rng.poisson(depth * rng.lognormal(-0.27**2 / 2, 0.27, len(oligos)))
    if label == "badread":
        reads, truth = nanopore.badread_reads(oligos, copies, out / "badread_tmp", seed=seed)
        (out / "badread_tmp" / "badread.fastq").unlink(missing_ok=True)
        return [s for _, s in reads], [t["oligo"] for t in truth]
    p = ChannelParams(p_sub=0.04, p_ins=0.025, p_del=0.025, mean_coverage=depth, coverage="lognormal",
                      lognormal_sigma=0.27, seed=seed)
    reads, _ = simulate(oligos, p, rng)
    return [r.seq for r in reads], [r.oligo for r in reads]


def cluster_metrics(reads: list[str], labels: list[int], primers, body_nt: int) -> dict:
    f, t = primers.forward, primers.tail
    bodies, lab = [], []
    for s, o in zip(reads, labels):
        x = locate_and_orient(s, f, t)
        if x is not None:
            bodies.append(x[0])
            lab.append(o)
    t0 = time.perf_counter()
    order = sorted(range(len(bodies)), key=lambda i: (len(bodies[i]) != body_nt, i))
    cl = kmer_cluster(bodies, order=order)
    t_cl = time.perf_counter() - t0
    multi = [m for m in cl if len(m) >= 2]
    impure = minority = 0
    per_oligo = defaultdict(list)
    for ci, m in enumerate(cl):
        c = Counter(lab[i] for i in m)
        if len(m) >= 2 and len(c) > 1:
            impure += 1
        minority += len(m) - c.most_common(1)[0][1]
        for o, n in c.items():
            per_oligo[o].append(n)
    split = sum(sum(n >= 2 for n in v) >= 2 for v in per_oligo.values())
    return {"reads_oriented": len(bodies), "clusters": len(cl), "clusters_ge2": len(multi),
            "merge_rate": impure / max(1, len(multi)), "minority_read_frac": minority / max(1, len(bodies)),
            "split_oligo_frac": split / max(1, len(per_oligo)),
            "clusters_per_oligo": float(np.mean([len(v) for v in per_oligo.values()])),
            "singleton_frac": sum(len(m) == 1 for m in cl) / max(1, len(cl)), "cluster_s": t_cl}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bytes", type=int, default=100_000)
    ap.add_argument("--depth", type=float, default=10.0)
    ap.add_argument("--bound", type=float, default=0.40)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--arms", default=",".join(a for a, _, _ in ARMS))
    ap.add_argument("--channels", default="badread,ids9")
    ap.add_argument("--out", default=str(ROOT / "results" / "phase3" / "cluster_scale"))
    args = ap.parse_args()
    out = Path(args.out)
    (out / "pools").mkdir(parents=True, exist_ok=True)
    primers = load_primers(ROOT / "configs" / "primers.yaml")
    data = random.Random(0).randbytes(args.bytes)
    rows = []
    for label, name, params in ARMS:
        if label not in args.arms.split(","):
            continue
        cfg = load_config(ROOT / "configs" / "default.yaml",
                          {"codec": {"name": name, "params": params}, "workers": args.workers})
        fa = out / "pools" / f"{label}.fa"
        t0 = time.perf_counter()
        res = encode_bytes(data, primers, cfg)
        enc_s = time.perf_counter() - t0
        write_fasta(res.oligos, str(fa))
        oligos = [s for _, s in res.oligos]
        lf, lt = len(primers.forward), len(primers.tail)
        body_nt = cfg["oligo_len"] - lf - lt
        t1 = time.perf_counter()
        md = min_distance([s[lf:-lt] for s in oligos], args.bound, args.workers)
        print(f"{label}: {len(oligos)} oligos, enc {enc_s:.0f}s, min body ed {md['min_ed']} "
              f"({md['min_ed_frac']}), pairs<=0.30L {md['pairs_le_0.3']}, <=0.36L {md['pairs_le_0.36']} "
              f"[{time.perf_counter() - t1:.0f}s]", flush=True)
        for ch in args.channels.split(","):
            reads, labels = channel_reads(ch, oligos, args.depth, 1, out)
            cm = cluster_metrics(reads, labels, primers, body_nt)
            dec = {}
            for st in ("erasure", "repair"):
                t2 = time.perf_counter()
                r = decode_reads(reads, primers, cfg["oligo_len"], workers=args.workers, strength=st)
                dec[st] = {"ok": r.ok, "multi_index": r.report.get("groups_multi_index"),
                           "valid": r.report.get("data_oligos_valid", r.report.get("droplets_valid")),
                           "repaired": r.report.get("repaired"), "error": r.report.get("error"),
                           "decode_s": time.perf_counter() - t2}
            row = {"arm": label, "channel": ch, "n_oligos": len(oligos), "encode_s": enc_s, **md, **cm,
                   **{f"{st}_{k}": v for st, d in dec.items() for k, v in d.items()}}
            rows.append(row)
            print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()
                              if k not in ("pairs",)}), flush=True)
            pd.DataFrame(rows).to_csv(out / "cluster_scale.csv", index=False)
    (out / "run.json").write_text(json.dumps(vars(args), indent=1))


if __name__ == "__main__":
    main()
