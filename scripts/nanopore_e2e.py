"""Phase 3 end-to-end check: encode -> read channel -> OUR decoder, per arm.

This is a pipeline validation on one payload, NOT the Phase 4 benchmark (no equal-nt
matching, no bisection, few subsample replicates).

Per arm:
1. Encode ``--bytes`` random bytes (seed 0) at default settings.
2. Draw lognormal copy numbers (sigma ``--sigma``, mean ``--max-depth``) per oligo and
   simulate with the selected Badread, uniform IDS, or squigulator -> Dorado channel.
3. **Grouping analysis** (DESIGN.md 0b.6) on all reads, using the truth table: per read,
   was it oriented, was the orientation right, and did its index field decode to the true
   index, to nothing (invalid repeat / wrong length), or to another index?
4. **Depth subsampling:** for each mean depth d, draw round(d * n_oligos) reads without
   replacement (``--reps`` replicates) and decode at both decoder strengths
   (``erasure``: CRC-discard + RS erasures; ``repair``: CRC-guided homopolymer repair + RS
   errors-and-erasures). This approximates a run of depth d from
   the same library (reads are a subsample of one simulated molecule set).
   Recorded: file recovery, data-oligo CRC failures, and post-consensus per-oligo error
   (edit distance of the decoded/attempted body vs the true body / body length; an index
   with no group scores 1.0 and is also counted as missing).

Channels: ``badread`` (read-level, calibrated to real R10.4.1 reads; no read-to-read
systematic errors), ``ids`` (uniform controlled channel), or ``squigulator`` (signal
simulation -> Dorado; miscalibrated, see docs/notes/channel_calibration.md).

Usage: .venv/bin/python scripts/nanopore_e2e.py --channel badread --out results/phase3/e2e_badread
"""

from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

import edlib
import numpy as np
import pandas as pd

from dnastore import nanopore
from dnastore.channel import ChannelParams, simulate
from dnastore.config import load_config
from dnastore.consensus import locate_and_orient
from dnastore.decoder import decode_reads
from dnastore.encoder import encode_bytes
from dnastore.metrics import pool_metrics
from dnastore.oligo import HEADER_INDEX_SPACE, INDEX_TRITS, decode_index
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


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(c - h), float(c + h)


def grouping(reads, truth: dict, oligos: list[tuple[int, str]], primers) -> tuple[dict, pd.DataFrame]:
    f, t = primers.forward, primers.tail
    valid = {idx for idx, _ in oligos}
    rows = []
    for rid, seq, _ in reads:
        oi, rev = truth[rid]
        true_idx = oligos[oi][0]
        o = locate_and_orient(seq, f, t)
        row = {"read_id": rid, "true_index": true_idx, "len": len(seq)}
        if o is None:
            row["outcome"] = "unoriented"
        else:
            body, is_rev, _ = o
            row["strand_ok"] = is_rev == rev
            row["body_len_err"] = len(body) - (len(oligos[oi][1]) - len(f) - len(t))
            idx = decode_index(body[:INDEX_TRITS], f[-1])
            if idx is None:
                row["outcome"] = "bad_index"
            elif idx == true_idx:
                row["outcome"] = "correct"
            elif idx in valid:
                row["outcome"] = "misassigned_valid"
            else:
                row["outcome"] = "misassigned_unused"
        rows.append(row)
    df = pd.DataFrame(rows)
    n = len(df)
    summ = {f"frac_{k}": v / n for k, v in df["outcome"].value_counts().items()}
    summ["n_reads"] = n
    summ["strand_ok_frac"] = float(df["strand_ok"].dropna().mean()) if "strand_ok" in df else float("nan")
    return summ, df


def per_oligo_error(res, oligos: list[tuple[int, str]], body_start: int, body_end: int) -> dict:
    """Post-consensus per-oligo error of the decoded/attempted body vs the true body.

    ``oligo_err_mean`` scores a missing oligo (no body) as 1.0; ``oligo_err_present_mean``
    averages only oligos that produced a body."""
    errs, present, missing, wrong = [], [], 0, 0
    for idx, seq in oligos:
        if idx < HEADER_INDEX_SPACE:
            continue
        true_body = seq[body_start:body_end]
        got = res.consensus.get(idx)
        if got is None:
            missing += 1
            errs.append(1.0)
            continue
        e = edlib.align(got, true_body, task="distance")["editDistance"] / len(true_body)
        errs.append(e)
        present.append(e)
        wrong += int(e > 0)
    e = np.array(errs)
    return {"oligo_err_mean": float(e.mean()), "oligo_err_nonzero_frac": float((e > 0).mean()),
            "oligo_err_present_mean": float(np.mean(present)) if present else float("nan"),
            "oligo_missing": missing, "oligo_wrong_body": wrong}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bytes", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-depth", type=float, default=40.0)
    ap.add_argument("--sigma", type=float, default=0.27, help="synthesis lognormal sigma (Gimpel 2023, Twist)")
    ap.add_argument("--depths", default="3,4,5,6,8,10,15,20,30")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--model", default="hac")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--arms", default=",".join(a for a, _, _ in ARMS))
    ap.add_argument("--strengths", default="erasure,repair")
    ap.add_argument("--channel", choices=["squigulator", "badread", "ids"], default="squigulator",
                    help="squigulator -> Dorado, calibrated Badread, or uniform IDS channel")
    ap.add_argument("--ids-rate", type=float, default=0.09,
                    help="total substitution + insertion + deletion rate for --channel ids; split equally")
    ap.add_argument("--out", default=str(ROOT / "results" / "phase3" / "e2e"))
    args = ap.parse_args()
    if args.channel == "squigulator" and not nanopore.tools_available(args.model):
        raise SystemExit("squigulator/Dorado not installed; run scripts/setup_nanopore_tools.sh")
    if not 0.0 <= args.ids_rate < 1.0:
        raise SystemExit("--ids-rate must be in [0, 1)")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    primers = load_primers(ROOT / "configs" / "primers.yaml")
    data = random.Random(args.seed).randbytes(args.bytes)
    depths = [float(d) for d in args.depths.split(",")]
    strengths = args.strengths.split(",")
    wanted = set(args.arms.split(","))
    grp_rows, dec_rows, enc_rows = [], [], []
    for label, name, params in ARMS:
        if label not in wanted:
            continue
        cfg = load_config(ROOT / "configs" / "default.yaml",
                          {"codec": {"name": name, "params": params}, "workers": args.workers})
        t0 = time.perf_counter()
        enc = encode_bytes(data, primers, cfg)
        enc_s = time.perf_counter() - t0
        pm = pool_metrics(enc)
        enc_rows.append({"arm": label, "encode_s": enc_s, **{k: pm[k] for k in
                         ("n_oligos", "n_header_oligos", "code_density", "effective_density", "total_nt")}})
        oligos = enc.oligos
        rng = np.random.default_rng([args.seed, len(enc_rows)])
        sig = args.sigma
        copies = rng.poisson(args.max_depth * rng.lognormal(-sig * sig / 2, sig, len(oligos)))
        run_dir = out / f"run_{label}"
        run_dir.mkdir(parents=True, exist_ok=True)
        if args.channel == "badread":
            t_ch = time.perf_counter()
            rl, tr = nanopore.badread_reads([s for _, s in oligos], copies, run_dir, seed=args.seed + 1)
            meta = {"t_dorado_s": 0.0, "t_badread_s": time.perf_counter() - t_ch}
            reads = [(rid, seq, "") for rid, seq in rl]
            truth = {t["read_id"]: (t["oligo"], t["reverse"]) for t in tr}
        elif args.channel == "ids":
            rate = args.ids_rate / 3.0
            params = ChannelParams(p_sub=rate, p_ins=rate, p_del=rate, mean_coverage=args.max_depth,
                                   coverage="lognormal", lognormal_sigma=sig, seed=args.seed + 1)
            t_ch = time.perf_counter()
            simulated, _ = simulate([s for _, s in oligos], params,
                                    np.random.default_rng(args.seed + 1))
            t_ids = time.perf_counter() - t_ch
            reads = [(f"ids_{i}", r.seq, "") for i, r in enumerate(simulated)]
            truth = {f"ids_{i}": (r.oligo, r.reverse) for i, r in enumerate(simulated)}
            meta = {"t_dorado_s": 0.0, "t_ids_s": t_ids, "n_reads": len(reads),
                    "ids_rate": args.ids_rate, "ids_sub": rate, "ids_ins": rate, "ids_del": rate}
        else:
            fq, meta = nanopore.run_pipeline([s for _, s in oligos], copies, run_dir, seed=args.seed + 1,
                                             model=args.model)
            reads = nanopore.read_fastq(fq)
            truth = {}
            for line in open(run_dir / "truth.tsv").readlines()[1:]:
                rid, oi, _, rev = line.split("\t")
                truth[rid] = (int(oi), rev.strip() == "1")
        gs, gdf = grouping(reads, truth, oligos, primers)
        gdf.to_csv(run_dir / "grouping.csv", index=False)
        grp_rows.append({"arm": label, **gs})
        print(f"{label:18s} oligos={len(oligos)} reads={len(reads)} enc={enc_s:.0f}s "
              f"channel={meta.get('t_badread_s', meta.get('t_ids_s', meta.get('t_dorado_s', 0))):.0f}s "
              f"grouping={ {k: round(v, 4) for k, v in gs.items() if k.startswith('frac')} }", flush=True)
        f, t = primers.forward, primers.tail
        seqs = [s for _, s, _ in reads]
        with open(run_dir / "oligos.fa", "w") as fh:
            for idx, s in oligos:
                fh.write(f">oligo_{idx}\n{s}\n")

        def row(d, r, st, res):
            rp = res.report
            ranks = rp.get("repair_hit_ranks") or []
            return {"arm": label, "depth": d, "rep": r, "strength": st, "ok": res.ok,
                    **per_oligo_error(res, oligos, len(f), cfg["oligo_len"] - len(t)),
                    **{k: rp.get(k) for k in ("data_crc_fail", "data_oligos_valid", "droplets_valid", "repaired",
                                              "repair_attempts", "repair_candidates_tried", "unverified_used",
                                              "rs_ee_columns", "rs_errors_mode", "rs_blocks_failed",
                                              "header_repaired", "error", "decode_s")},
                    "repair_rank_max": max(ranks, default=None),
                    "repair_edits_2": sum(e == 2 for e in rp.get("repair_hit_edits") or [])}

        for st in strengths:
            dec_rows.append(row(len(seqs) / len(oligos), -1, st,
                                decode_reads(seqs, primers, cfg["oligo_len"], workers=args.workers, strength=st)))
        for d in depths:
            n_take = int(round(d * len(oligos)))
            if n_take > len(seqs):
                continue
            for r in range(args.reps):
                sub_rng = np.random.default_rng([args.seed, int(d * 100), r])
                pick = sub_rng.choice(len(seqs), n_take, replace=False)
                sub_reads = [seqs[i] for i in pick]
                for st in strengths:
                    res = decode_reads(sub_reads, primers, cfg["oligo_len"], workers=args.workers, strength=st)
                    dec_rows.append(row(d, r, st, res))
        sub = pd.DataFrame([x for x in dec_rows if x["arm"] == label])
        print(sub.groupby(["strength", "depth"]).agg(ok=("ok", "mean"), err=("oligo_err_mean", "mean"),
                                                     wrong=("oligo_wrong_body", "mean"),
                                                     repaired=("repaired", "mean")).round(4).to_string(), flush=True)
    pd.DataFrame(enc_rows).to_csv(out / "encode.csv", index=False)
    pd.DataFrame(grp_rows).to_csv(out / "grouping_summary.csv", index=False)
    dec = pd.DataFrame(dec_rows)
    dec.to_csv(out / "decode.csv", index=False)
    rows = []
    for (arm, st, d), g in dec[dec.rep >= 0].groupby(["arm", "strength", "depth"]):
        lo, hi = wilson(int(g.ok.sum()), len(g))
        rows.append({"arm": arm, "strength": st, "depth": d, "n": len(g), "recovered": int(g.ok.sum()), "wilson_lo": lo,
                     "wilson_hi": hi, "oligo_err_mean": g.oligo_err_mean.mean(),
                     "oligo_missing_mean": g.oligo_missing.mean(), "data_crc_fail_mean": g.data_crc_fail.mean()})
    pd.DataFrame(rows).to_csv(out / "decode_summary.csv", index=False)
    (out / "run.json").write_text(json.dumps({**vars(args), "versions": nanopore.tool_versions()}, indent=2,
                                             default=str))


if __name__ == "__main__":
    main()
