"""Validate cheap structure proxies against primer-site accessibility (primary metric).

Builds >= 1000 oligos with the real primers:
  * ``quat``: index field (15 rotating-ternary nt) + 145 nt i.i.d. quaternary (whitened payload),
  * ``tern``: 160 nt rotating ternary (goldman / alphabet-B-like body),
  * ``adv``:  quaternary body with an inserted 6-14 nt substring of R or rc(F) (complementary to
              a primer site) to span the low-accessibility range.

For each oligo it computes accessibility (full partition function, both strands) and MFE at
60 and 37 C, limited-span MFE (reference only), GC fraction, and the incremental proxy's
components for several k sets. Weights for combined proxies are chosen on the EVEN half
of the oligos; all reported rho values come from the ODD (held-out) half. Spearman rho uses
average ranks and a bootstrap (2000 resamples) for 95% CIs.

Usage: .venv/bin/python scripts/validate_structure_proxy.py [--n-quat 450 --n-tern 450 --n-adv 110 --workers 8]
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from dnastore.dnautil import BASES, revcomp, rotate_encode
from dnastore.primers import load_primers
from dnastore.structure import (
    StructureConfig,
    accessibility,
    limited_span_mfe,
    local_access,
    mfe_at,
    proxy_components,
)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "phase2" / "structure"
KSETS = {"k4-6": (4, 5, 6), "k5-6": (5, 6), "k6-7": (6, 7), "k4-7": (4, 5, 6, 7)}
WINDOWS = (40, 60, 80, 100)
SALT = StructureConfig().salt_molar


def build_oligos(n_quat: int, n_tern: int, n_adv: int, seed: int) -> list[dict]:
    P = load_primers(ROOT / "configs" / "primers.yaml")
    F, tail, R = P.forward, P.tail, P.reverse
    rng = np.random.default_rng(seed)
    out = []

    def quat_body():
        idx = rotate_encode(list(rng.integers(0, 3, 15)), F[-1])
        return idx + "".join(BASES[b] for b in rng.integers(0, 4, 145))

    for _ in range(n_quat):
        out.append({"type": "quat", "oligo": F + quat_body() + tail})
    for _ in range(n_tern):
        out.append({"type": "tern", "oligo": F + rotate_encode(list(rng.integers(0, 3, 160)), F[-1]) + tail})
    rcF = revcomp(F)
    for _ in range(n_adv):
        body = quat_body()
        src = R if rng.random() < 0.5 else rcF
        ln = int(rng.integers(6, 15))
        s = int(rng.integers(0, len(src) - ln + 1))
        frag = src[s : s + ln]
        pos = int(rng.integers(15, len(body) - ln))
        body = body[:pos] + frag + body[pos + ln :]
        out.append({"type": "adv", "oligo": F + body + tail})
    return out


def _metrics(job: tuple[int, dict, int, int]) -> dict:
    i, rec, f_len, r_len = job
    s = rec["oligo"]
    row = {"i": i, "type": rec["type"], "oligo": s}
    tm = {}
    for T in (60.0, 37.0):
        tag = f"_T{int(T)}"
        t0 = time.perf_counter()
        for k, v in accessibility(s, f_len, r_len, T, SALT).items():
            row[k + tag] = v
        tm["access" + tag] = time.perf_counter() - t0
        t0 = time.perf_counter()
        row["mfe_min" + tag] = min(mfe_at(s, T, SALT), mfe_at(revcomp(s), T, SALT))
        tm["mfe_both" + tag] = time.perf_counter() - t0
    for span in (30, 50):
        t0 = time.perf_counter()
        row[f"lsmfe{span}_T60"] = min(limited_span_mfe(s, 60.0, SALT, span), limited_span_mfe(revcomp(s), 60.0, SALT, span))
        tm[f"lsmfe{span}"] = time.perf_counter() - t0
    row["gc"] = (s.count("G") + s.count("C")) / len(s)
    for name, ks in KSETS.items():
        t0 = time.perf_counter()
        c = proxy_components(s, f_len, r_len, ks)
        tm[f"proxy_{name}"] = (time.perf_counter() - t0) / 2  # components computes the proxy twice
        row[f"self_{name}"], row[f"site_{name}"] = c["self"], c["site"]
    row.update({f"time_{k}": v for k, v in tm.items()})
    return row


def _local(job: tuple[str, int, int]) -> dict:
    s, f_len, r_len = job
    row = {}
    for w in WINDOWS:
        t0 = time.perf_counter()
        la = local_access(s, f_len, r_len, w, 60.0, SALT)
        row[f"time_local{w}"] = time.perf_counter() - t0
        row[f"local{w}_min_anchor_T60"] = la["local_min_anchor"]
        row[f"local{w}_min_site_T60"] = la["local_min_site"]
    return row


def rankdata(x: np.ndarray) -> np.ndarray:
    return pd.Series(x).rank(method="average").to_numpy()


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        boots.append(spearman(x[idx], y[idx]))
    boots = np.array(boots)
    boots = boots[~np.isnan(boots)]
    return spearman(x, y), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-quat", type=int, default=450)
    ap.add_argument("--n-tern", type=int, default=450)
    ap.add_argument("--n-adv", type=int, default=110)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20261001)
    ap.add_argument("--analyze-only", action="store_true", help="Reuse proxy_validation.csv (skip the expensive pass).")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    P = load_primers(ROOT / "configs" / "primers.yaml")
    f_len, r_len = len(P.forward), len(P.reverse)
    csv = OUT / "proxy_validation.csv"
    wall = None
    if args.analyze_only:
        df = pd.read_csv(csv)
    else:
        recs = build_oligos(args.n_quat, args.n_tern, args.n_adv, args.seed)
        jobs = [(i, r, f_len, r_len) for i, r in enumerate(recs)]
        t0 = time.perf_counter()
        with ProcessPoolExecutor(args.workers) as ex:
            rows = list(ex.map(_metrics, jobs, chunksize=4))
        wall = time.perf_counter() - t0
        df = pd.DataFrame(rows).sort_values("i")
    if f"local{WINDOWS[-1]}_min_anchor_T60" not in df.columns:
        with ProcessPoolExecutor(args.workers) as ex:
            loc = list(ex.map(_local, [(s, f_len, r_len) for s in df.oligo], chunksize=8))
        df = pd.concat([df.reset_index(drop=True), pd.DataFrame(loc)], axis=1)
    df.to_csv(csv, index=False)

    train, test = df[df.i % 2 == 0], df[df.i % 2 == 1]
    targets = ["min_site_access_T60", "min_anchor_access_T60", "min_site_access_T37", "min_anchor_access_T37"]

    # Combined proxies: score = self + lam * site; choose lam on train vs primary target.
    combos = {}
    for name in KSETS:
        best = None
        for lam in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 1e6):
            sc = train[f"self_{name}"] + lam * train[f"site_{name}"]
            r = spearman(sc.to_numpy(), train["min_anchor_access_T60"].to_numpy())
            if best is None or r < best[1]:  # most negative = best (more structure -> less access)
                best = (lam, r)
        combos[name] = best[0]
        df[f"combo_{name}"] = df[f"self_{name}"] + best[0] * df[f"site_{name}"]
    test = df[df.i % 2 == 1]

    candidates = [c for c in df.columns if c.startswith(("self_k", "site_k", "combo_k"))] + [
        f"local{w}_min_anchor_T60" for w in WINDOWS] + [
        "lsmfe30_T60", "lsmfe50_T60", "mfe_min_T60", "mfe_min_T37", "gc"]
    res = {}
    for c in candidates:
        res[c] = {}
        for t in targets:
            rho, lo, hi = spearman_ci(test[c].to_numpy(float), test[t].to_numpy(float))
            res[c][t] = {"rho": rho, "ci95": [lo, hi]}
        # per-type on the primary target (held-out)
        for typ in ("quat", "tern", "adv"):
            sub = test[test.type == typ]
            res[c][f"min_anchor_access_T60|{typ}"] = spearman(sub[c].to_numpy(float), sub["min_anchor_access_T60"].to_numpy(float))

    timing = {c[5:]: {"mean_ms": 1000 * df[c].mean()} for c in df.columns if c.startswith("time_")}

    # Threshold calibration from accessibility (DESIGN Part 0 item 2; see docs/notes/structure.md).
    # Pessimistic model: per-cycle efficiency proportional to equilibrium 3'-anchor accessibility a.
    # Relative yield after n cycles vs a fully accessible oligo (efficiency 1): ((1+a)/2)**n.
    # Require <= 2-fold under-representation after n = 25 cycles: a >= 2*0.5**(1/25) - 1.
    # (i) absolute: reference = fully accessible site (a_ref = 1).
    # (ii) relative: reference = median min-anchor accessibility of UNCONSTRAINED random bodies
    #      (quat + tern), i.e. "no more than `fold` x under-represented vs a typical oligo".
    #      a_min = (1 + a_ref) * fold**(-1/n) - 1.
    unc = df[df.type.isin(["quat", "tern"])]["min_anchor_access_T60"]
    a_ref_rel = float(unc.median())
    thr = {}
    for n in (20, 25, 30):
        for fold in (2.0, 4.0):
            thr[f"abs_n{n}_fold{int(fold)}"] = 2 * (1 / fold) ** (1 / n) - 1
            thr[f"rel_n{n}_fold{int(fold)}"] = (1 + a_ref_rel) * (1 / fold) ** (1 / n) - 1
    a_star = thr["rel_n25_fold2"]
    dist = {}
    for typ in ("quat", "tern", "adv", "all"):
        sub = df if typ == "all" else df[df.type == typ]
        d = {}
        for t in ("min_anchor_access_T60", "min_site_access_T60", "min_anchor_access_T37"):
            v = sub[t].to_numpy()
            d[t] = {"p1": float(np.percentile(v, 1)), "p5": float(np.percentile(v, 5)), "p25": float(np.percentile(v, 25)),
                    "median": float(np.median(v)), "mean": float(v.mean())}
        for name, a in thr.items():
            d[f"frac_fail_anchor_T60_{name}"] = float((sub["min_anchor_access_T60"] < a).mean())
        d["n"] = int(len(sub))
        dist[typ] = d

    summary = {
        "n_oligos": len(df), "n_train": len(train), "n_test": len(test), "workers": args.workers,
        "wall_s": wall, "seed": args.seed, "salt_molar": SALT,
        "combo_lambda_selected_on_train": combos,
        "spearman_heldout": res, "timing_ms_per_oligo": timing,
        "threshold_candidates": thr, "a_ref_relative_median_unconstrained": a_ref_rel,
        "proposed_threshold_min_anchor_access_T60": a_star, "proposed_rule": "rel_n25_fold2",
        "distributions": dist,
    }
    (OUT / "proxy_validation_summary.json").write_text(json.dumps(summary, indent=2))
    # concise console table
    print(f"n={len(df)} wall={wall if wall is None else round(wall)}s")
    for c in candidates:
        r = res[c]["min_anchor_access_T60"]
        s = res[c]["min_site_access_T60"]
        print(f"{c:22s} anchor60 rho={r['rho']:+.3f} [{r['ci95'][0]:+.3f},{r['ci95'][1]:+.3f}]  "
              f"site60 rho={s['rho']:+.3f} [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]  "
              f"quat={res[c]['min_anchor_access_T60|quat']:+.2f} tern={res[c]['min_anchor_access_T60|tern']:+.2f}")
    print("lambda:", combos)
    print("timing:", {k: round(v["mean_ms"], 2) for k, v in timing.items()})
    print("threshold a* =", a_star)
    for typ, d in dist.items():
        print(typ, {k: round(v, 3) for k, v in d["min_anchor_access_T60"].items()},
              {k[len("frac_fail_anchor_T60_"):]: round(v, 3) for k, v in d.items() if k.startswith("frac_fail")})


if __name__ == "__main__":
    main()
