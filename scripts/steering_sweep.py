"""Phase 2 characterisation: what each arm's encoder produces (NOT the Phase 4 benchmark).

For one fixed random payload, encode with every arm and P value, then measure per data
oligo: hard-constraint violations, rerolls, GC, max run, confusability (R10 delta=0.05:
mean/max over both strands), Whritenour adjacent contrast (R9 pA, tau=4), the stem proxy,
and structure (primer-site/3'-anchor accessibility and MFE at 60 and 37 deg C).
No channel and no decoding under noise: those come in Phases 3-4.

Usage: .venv/bin/python scripts/steering_sweep.py --bytes 4000 --workers 8 --out results/phase2/sweep
"""

from __future__ import annotations

import argparse
import json
import random
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from dnastore.config import load_config
from dnastore.encoder import encode_bytes
from dnastore.metrics import pool_metrics
from dnastore.oligo import HEADER_INDEX_SPACE
from dnastore.primers import load_primers
from dnastore.screening import ScreeningConfig, screen_oligo

ROOT = Path(__file__).resolve().parent.parent
ARMS = [("naive2bit", {}), ("naive2bit", {"whiten": True}), ("goldman", {})]
for _alpha in ("A", "C", "B"):
    for _P in (4, 6, 8, 12, 16, 0):
        ARMS.append(("steering", {"P": _P, "alphabet": _alpha}))

_CTX: dict = {}


def _init() -> None:
    from dnastore import kmer

    _CTX["r10"] = kmer.load_model("r10.4.1_9mer")
    _CTX["r9"] = kmer.load_model("r9.4.1_6mer")
    _CTX["conf"] = kmer.confusability_table(_CTX["r10"], delta=0.05)
    _CTX["conf9"] = kmer.confusability_table(_CTX["r9"], c=1.0)
    _CTX["primers"] = load_primers(ROOT / "configs" / "primers.yaml")


def _measure(args: tuple) -> dict:
    from dnastore import kmer, structure

    if not _CTX:
        _init()
    oligo, scfg_d = args
    p = _CTX["primers"]
    lf, lr = len(p.forward), len(p.reverse)
    row = {}
    scr = screen_oligo(oligo, p, ScreeningConfig.from_dict(scfg_d), compute_mfe=False)
    row.update(gc=scr["gc"], max_hp=scr["max_hp"], gc_win_min=scr["gc_win_min"], gc_win_max=scr["gc_win_max"],
               primer_hits=scr["primer_hits"], anchor_hits=scr["primer_anchor_hits"],
               violations="|".join(scr["violations"]))
    c = kmer.oligo_confusability(oligo, _CTX["conf"], _CTX["r10"].k)
    floor = float(np.quantile(_CTX["conf"], 0.75))
    codes = kmer.encode_kmers(oligo, _CTX["r10"].k)
    rc = kmer.rc_code_table(_CTX["r10"].k)[codes]
    tail = float(np.maximum(_CTX["conf"][codes] - floor, 0).sum() + np.maximum(_CTX["conf"][rc] - floor, 0).sum())
    c9 = kmer.oligo_confusability(oligo, _CTX["conf9"], _CTX["r9"].k)
    row.update(conf_mean=c["mean"], conf_max=c["max"], conf_tail=tail, conf9_mean=c9["mean"], conf9_max=c9["max"])
    adj = kmer.adjacent_contrast(oligo, _CTX["r9"], tau=4.0)
    row.update(adj_r9_min=adj["min"], adj_r9_frac_below=adj["frac_below_tau"])
    row["stem_proxy"] = structure.proxy_score(oligo, lf, lr, structure.ProxyParams((4, 5, 6), 3, 1.0, 8.0))
    row["local_min_anchor_W80"] = structure.local_access(oligo, lf, lr, 80, 60.0, 0.05)["local_min_anchor"]
    row.update(structure.structure_report(oligo, lf, lr))
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bytes", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=str(ROOT / "results" / "phase2" / "sweep"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    primers = load_primers(ROOT / "configs" / "primers.yaml")
    data = random.Random(args.seed).randbytes(args.bytes)
    rows, summaries = [], []
    with ProcessPoolExecutor(args.workers, initializer=_init) as ex:
        for name, params in ARMS:
            label = name + ("-whiten" if params.get("whiten") else "") + (
                f"-{params['alphabet']}-P{params['P'] or 'inf'}" if name == "steering" else "")
            cfg = load_config(ROOT / "configs" / "default.yaml",
                              {"codec": {"name": name, "params": params}, "workers": args.workers})
            t0 = time.perf_counter()
            res = encode_bytes(data, primers, cfg)
            enc_s = time.perf_counter() - t0
            data_oligos = [(i, s) for i, s in res.oligos if i >= HEADER_INDEX_SPACE]
            info = {d["index"]: d for d in res.payload_info}
            t1 = time.perf_counter()
            meas = list(ex.map(_measure, [(s, cfg["screening"]) for _, s in data_oligos], chunksize=4))
            for (idx, _), m in zip(data_oligos, meas):
                rows.append({"arm": label, "index": idx, "seed": info[idx]["seed"], "tries": info[idx]["tries"],
                             "recorded_violations": "|".join(info[idx]["violations"]), **m})
            pm = pool_metrics(res)
            summaries.append({"arm": label, "encode_s": enc_s, "measure_s": time.perf_counter() - t1, **pm})
            print(f"{label:22s} oligos={pm['n_oligos']:4d} code_density={pm['code_density']:.3f} "
                  f"enc={enc_s:.1f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "per_oligo.csv", index=False)
    agg = df.groupby("arm").agg(
        frac_violating=("violations", lambda v: float((v != "").mean())),
        max_tries=("tries", "max"),
        max_hp=("max_hp", "max"), conf_tail=("conf_tail", "mean"), conf_max_mean=("conf_max", "mean"),
        conf9_max_mean=("conf9_max", "mean"),
        adj_r9_frac_below=("adj_r9_frac_below", "mean"), stem_proxy=("stem_proxy", "median"),
        min_anchor_access_T60=("min_anchor_access_T60", "median"),
        frac_access_686=("min_anchor_access_T60", lambda x: float((x >= 0.686).mean())),
        min_site_access_T60=("min_site_access_T60", "median"),
        mfe_rc_T60=("mfe_rc_T60", "median"), mfe_fwd_T37=("mfe_fwd_T37", "median"),
    ).reset_index()
    summ = pd.DataFrame(summaries)[["arm", "n_oligos", "code_density", "effective_density", "seed_exhaustion_rate",
                                    "mean_tries", "encode_s"]]
    table = summ.merge(agg, on="arm")
    table.to_csv(out / "summary.csv", index=False)
    (out / "run.json").write_text(json.dumps({"bytes": args.bytes, "seed": args.seed, "arms": [a for a, _ in ARMS]}))
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
