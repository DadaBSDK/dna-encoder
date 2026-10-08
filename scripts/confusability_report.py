"""Confusability report (Phase 2): per-k-mer distributions and per-oligo statistics.

Writes to results/phase2/confusability/:
  kmer_table_stats.csv    distribution of conf(x) over all k-mers, per model/parameter
  oligo_metrics.csv       per-oligo metrics for random whitened quaternary vs rotating-ternary bodies
  summary.json            group summaries and Spearman correlations

Usage: python scripts/confusability_report.py [--n 1000] [--seed 0]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dnastore.dnautil import rotate_encode
from dnastore.kmer import adjacent_contrast, confusability_table, load_model, oligo_confusability
from dnastore.primers import load_primers

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "phase2" / "confusability"
R10_DELTAS = [0.01, 0.025, 0.05, 0.1, 0.2]
R9_CS = [0.5, 1.0, 2.0]
BODY_NT = 160


def spearman(a: pd.Series, b: pd.Series) -> float:
    """Spearman rho as Pearson correlation of average ranks (no scipy dependency)."""
    return float(a.rank().corr(b.rank()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    r10, r9 = load_model("r10.4.1_9mer"), load_model("r9.4.1_6mer")
    tables = {f"r10_d{d}": (confusability_table(r10, delta=d), 9) for d in R10_DELTAS}
    tables.update({f"r9_c{c}": (confusability_table(r9, c=c), 6) for c in R9_CS})

    rows = []
    for name, (t, _) in tables.items():
        q = np.percentile(t, [5, 25, 50, 75, 95])
        rows.append({"table": name, "mean": t.mean(), "std": t.std(), "frac_zero": (t == 0).mean(),
                     "p5": q[0], "p25": q[1], "p50": q[2], "p75": q[3], "p95": q[4], "max": t.max()})
    pd.DataFrame(rows).to_csv(OUT / "kmer_table_stats.csv", index=False)

    primers = load_primers(ROOT / "configs" / "primers.yaml")
    rng = np.random.default_rng(args.seed)
    recs = []
    for kind in ("whitened_quaternary", "rotating_ternary"):
        for i in range(args.n):
            if kind == "whitened_quaternary":
                body = "".join("ACGT"[b] for b in rng.integers(0, 4, BODY_NT))
            else:
                body = rotate_encode(list(rng.integers(0, 3, BODY_NT)), primers.forward[-1])
            oligo = primers.forward + body + primers.tail
            r = {"kind": kind, "i": i}
            for name, (t, k) in tables.items():
                oc = oligo_confusability(oligo, t, k)
                r[f"{name}_mean"], r[f"{name}_max"] = oc["mean"], oc["max"]
            a9 = adjacent_contrast(oligo, r9, tau=4.0)
            a10 = adjacent_contrast(oligo, r10, tau=0.5)
            r.update(adj_r9_min=a9["min"], adj_r9_mean=a9["mean"], adj_r9_frac_le4pA=a9["frac_below_tau"],
                     adj_r10_min=a10["min"], adj_r10_mean=a10["mean"], adj_r10_frac_le0p5=a10["frac_below_tau"])
            recs.append(r)
    df = pd.DataFrame(recs)
    df.to_csv(OUT / "oligo_metrics.csv", index=False)

    metric_cols = [c for c in df.columns if c not in ("kind", "i")]
    summary = {
        "n_per_kind": args.n, "seed": args.seed, "body_nt": BODY_NT,
        "models": {m.name: {"sha256": m.sha256, "units": m.units, "k": m.k} for m in (r10, r9)},
        "group_means": df.groupby("kind")[metric_cols].mean().to_dict(),
        "spearman_all_oligos": {},
    }
    pairs = [
        ("r10_d0.05_mean", "r9_c1.0_mean"),
        ("r10_d0.05_mean", "adj_r9_frac_le4pA"),
        ("r10_d0.05_mean", "adj_r10_frac_le0p5"),
        ("r10_d0.05_mean", "adj_r10_mean"),
        ("r9_c1.0_mean", "adj_r9_frac_le4pA"),
        ("r9_c1.0_mean", "adj_r9_mean"),
        ("r10_d0.05_mean", "r10_d0.1_mean"),
        ("r10_d0.01_mean", "r10_d0.2_mean"),
    ]
    for a, b in pairs:
        summary["spearman_all_oligos"][f"{a}~{b}"] = spearman(df[a], df[b])
        for kind, g in df.groupby("kind"):
            summary.setdefault("spearman_within_kind", {}).setdefault(kind, {})[f"{a}~{b}"] = spearman(g[a], g[b])
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(pd.DataFrame(rows).to_string(index=False))
    print(df.groupby("kind")[["r10_d0.05_mean", "r10_d0.05_max", "r9_c1.0_mean", "adj_r9_mean",
                              "adj_r9_frac_le4pA", "adj_r10_mean"]].mean().to_string())
    print(json.dumps(summary["spearman_all_oligos"], indent=1))
    print(json.dumps(summary["spearman_within_kind"], indent=1))


if __name__ == "__main__":
    main()
