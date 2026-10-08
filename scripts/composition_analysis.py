"""Why does rotating ternary fold more stably than whitened quaternary?

Tests the candidate result "homopolymer constraints carry a hidden structure cost"
(DESIGN.md Part 0.9). The analysis is seeded and writes to
``results/phase2/composition/``.

Sources (all first-order Markov chains with a uniform stationary distribution, 160 nt):
  rot      rotating ternary with uniform trits (dnastore.dnautil.rotate_encode); P(repeat) = 0
  p0.05 .. p0.15   P(next == prev) = r, otherwise uniform over the other 3 bases
  iid      i.i.d. uniform quaternary (= r = 0.25, the "whitened" payload model)
  iid_hp3  i.i.d. quaternary rejection-filtered to max homopolymer <= 3 (the realistic constraint)
  shiftban control: forbids a -> next(a) (A->C, C->G, G->T, T->A) instead of a -> a.
           It has the same transition collision probability (sum_b T(a,b)^2 = 1/3) as rot,
           and is reverse-complement symmetric, but it ALLOWS homopolymers. If it folds like
           rot, the cause is the reduced per-base entropy, not the absence of repeats.

Theory: for a reverse-complement-symmetric chain, P(rc(w)) = P(w), so the probability
that two distant k-windows are mutually reverse-complementary is
sum_w P(w)^2 = 1/4 * c^(k-1), where c = sum_b T(a,b)^2 (1/4 for iid, 1/3 for rot).
Stacking: a helix formed by w and rc(w) has dinucleotide steps distributed proportionally
to T(a,b)^2, i.e. uniform over the allowed transitions. Its expected nearest-neighbour
dG37 per step is compared using SantaLucia (1998) unified parameters
(PNAS 95:1460-1465, doi:10.1073/pnas.95.4.1460).
"""

from __future__ import annotations

import json
import random
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from dnastore.dnautil import BASES, revcomp, rotate_encode  # noqa: E402

OUT = ROOT / "results" / "phase2" / "composition"
L = 160
N = 2000
N_PRIMER = 300
SEED = 20261001
TEMPS = (37.0, 60.0)
KS = (4, 5, 6, 7, 8)
MIN_LOOP = 3

# SantaLucia 1998 unified NN dG37 (kcal/mol), keyed by the 5'->3' top-strand dinucleotide.
NN_DG37 = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}
DINUCS = ["".join(p) for p in product(BASES, BASES)]
SELF_COMP = {"AT", "TA", "CG", "GC"}
NEXT = {b: BASES[(i + 1) % 4] for i, b in enumerate(BASES)}


# --------------------------------------------------------------------------- sources


def transition(kind: str, r: float = 0.25) -> dict[str, dict[str, float]]:
    if kind == "repeat":
        return {a: {b: (r if a == b else (1 - r) / 3) for b in BASES} for a in BASES}
    if kind == "shiftban":
        return {a: {b: (0.0 if b == NEXT[a] else 1 / 3) for b in BASES} for a in BASES}
    raise ValueError(kind)


def sample_markov(rng: random.Random, T: dict, n: int, first: str | None = None) -> str:
    s = [first or rng.choice(BASES)]
    for _ in range(n - 1):
        row = T[s[-1]]
        s.append(rng.choices(BASES, weights=[row[b] for b in BASES])[0])
    return "".join(s)


def max_run(s: str) -> int:
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def make_sources(rng: random.Random) -> dict[str, list[str]]:
    src: dict[str, list[str]] = {}
    src["rot"] = [rotate_encode([rng.randrange(3) for _ in range(L)], rng.choice(BASES)) for _ in range(N)]
    for r in (0.05, 0.10, 0.15):
        T = transition("repeat", r)
        src[f"p{r:.2f}"] = [sample_markov(rng, T, L) for _ in range(N)]
    src["iid"] = ["".join(rng.choice(BASES) for _ in range(L)) for _ in range(N)]
    hp3 = []
    while len(hp3) < N:
        s = "".join(rng.choice(BASES) for _ in range(L))
        if max_run(s) <= 3:
            hp3.append(s)
    src["iid_hp3"] = hp3
    T = transition("shiftban")
    src["shiftban"] = [sample_markov(rng, T, L) for _ in range(N)]
    return src


P_REPEAT = {"rot": 0.0, "p0.05": 0.05, "p0.10": 0.10, "p0.15": 0.15, "iid": 0.25}


# --------------------------------------------------------------------------- theory


def theory(T: dict) -> dict:
    c = float(np.mean([sum(T[a][b] ** 2 for b in BASES) for a in BASES]))  # same for every a here
    di = {a + b: 0.25 * T[a][b] for a in BASES for b in BASES}
    w = {a + b: T[a][b] ** 2 for a in BASES for b in BASES}
    z = sum(w.values())
    helix_dg = sum(NN_DG37[d] * w[d] / z for d in DINUCS)
    out = {
        "collision_c": c,
        "selfcomp_step_frac": sum(di[d] for d in SELF_COMP),
        "strong_step_frac_CG_GC_GG_CC": sum(di[d] for d in ("CG", "GC", "GG", "CC")),
        "mean_dG_per_step_sequence": sum(NN_DG37[d] * di[d] for d in DINUCS),
        "mean_dG_per_step_in_matched_helix": helix_dg,
        "dinuc": di,
    }
    for k in KS:
        p_pair = 0.25 * c ** (k - 1)
        nw = L - k + 1
        n_pairs = sum(max(0, nw - (i + k + MIN_LOOP)) for i in range(nw))
        out[f"exp_comp_pairs_k{k}"] = n_pairs * p_pair
    return out


# --------------------------------------------------------------------------- empirical composition


def comp_pairs(s: str, k: int) -> int:
    pos: dict[str, list[int]] = {}
    for i in range(len(s) - k + 1):
        pos.setdefault(s[i : i + k], []).append(i)
    n = 0
    for i in range(len(s) - k + 1):
        for j in pos.get(revcomp(s[i : i + k]), ()):
            if j >= i + k + MIN_LOOP:
                n += 1
    return n


def composition(seqs: list[str]) -> dict:
    c = Counter()
    for s in seqs:
        c.update(s[i : i + 2] for i in range(len(s) - 1))
    tot = sum(c.values())
    di = {d: c[d] / tot for d in DINUCS}
    out = {
        "gc": float(np.mean([(s.count("G") + s.count("C")) / len(s) for s in seqs])),
        "max_run_mean": float(np.mean([max_run(s) for s in seqs])),
        "selfcomp_step_frac": sum(di[d] for d in SELF_COMP),
        "strong_step_frac_CG_GC_GG_CC": sum(di[d] for d in ("CG", "GC", "GG", "CC")),
        "mean_dG_per_step_sequence": sum(NN_DG37[d] * di[d] for d in DINUCS),
        "dinuc": di,
    }
    for k in KS:
        out[f"obs_comp_pairs_k{k}"] = float(np.mean([comp_pairs(s, k) for s in seqs]))
    return out


# --------------------------------------------------------------------------- folding


def _fold(args: tuple[str, float]) -> tuple[float, str]:
    import RNA

    if not getattr(_fold, "_loaded", False):
        RNA.params_load_DNA_Mathews2004()
        _fold._loaded = True
    seq, temp = args
    md = RNA.md()
    md.temperature = temp
    st, e = RNA.fold_compound(seq, md).mfe()
    return float(e), st


def structure_stats(seq: str, st: str) -> dict:
    stack, pair = [], {}
    for i, ch in enumerate(st):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            pair[j] = i
    steps = []
    for i, j in pair.items():
        if pair.get(i + 1) == j - 1:  # stacked pair (i,j),(i+1,j-1)
            steps.append(seq[i : i + 2])
    paired = 2 * len(pair) / len(seq)
    dg = [NN_DG37[d] for d in steps]
    return {
        "paired_frac": paired,
        "n_stacks": len(steps),
        "stack_dG_mean": float(np.mean(dg)) if dg else np.nan,
        "gc_stack_frac": float(np.mean([d in ("CG", "GC", "GG", "CC") for d in steps])) if steps else np.nan,
    }


def boot_ci(x: np.ndarray, stat=np.median, n=2000, seed=0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    v = [stat(rng.choice(x, size=len(x), replace=True)) for _ in range(n)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def main() -> None:
    import pandas as pd

    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    src = make_sources(rng)

    # theory
    th = {"rot": theory(transition("repeat", 0.0)), "iid": theory(transition("repeat", 0.25)),
          "shiftban": theory(transition("shiftban"))}
    for r in (0.05, 0.10, 0.15):
        th[f"p{r:.2f}"] = theory(transition("repeat", r))

    comp = {name: composition(seqs) for name, seqs in src.items()}

    # folding
    jobs = [(name, i, s, t) for t in TEMPS for name, seqs in src.items() for i, s in enumerate(seqs)]
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(_fold, [(s, t) for _, _, s, t in jobs], chunksize=64))
    rows = []
    for (name, i, s, t), (e, st) in zip(jobs, res):
        rows.append({"source": name, "i": i, "temp": t, "mfe": e, **structure_stats(s, st),
                     **{f"cp_k{k}": None for k in ()}})
    df = pd.DataFrame(rows)
    # per-sequence complementary 6-mer pairs (for a within-source regression)
    cp6 = {(name, i): comp_pairs(s, 6) for name, seqs in src.items() for i, s in enumerate(seqs)}
    df["cp_k6"] = [cp6[(n, i)] for n, i in zip(df.source, df.i)]
    df.to_csv(OUT / "mfe_per_sequence.csv", index=False)

    # with primers (37 C only)
    from dnastore.primers import load_primers

    P = load_primers(ROOT / "configs" / "primers.yaml")
    prng = random.Random(SEED + 1)
    pseqs = {
        "iid": [P.forward + "".join(prng.choice(BASES) for _ in range(L)) + P.tail for _ in range(N_PRIMER)],
        "rot": [P.forward + rotate_encode([prng.randrange(3) for _ in range(L)], P.forward[-1]) + P.tail for _ in range(N_PRIMER)],
    }
    with ProcessPoolExecutor(8) as ex:
        pres = {k: [e for e, _ in ex.map(_fold, [(s, 37.0) for s in v], chunksize=16)] for k, v in pseqs.items()}

    # summaries
    summ = []
    for (name, t), g in df.groupby(["source", "temp"]):
        m = g.mfe.to_numpy()
        lo, hi = boot_ci(m)
        summ.append({
            "source": name, "temp": t, "p_repeat": P_REPEAT.get(name), "n": len(m),
            "mfe_median": float(np.median(m)), "mfe_median_ci_lo": lo, "mfe_median_ci_hi": hi,
            "mfe_mean": float(m.mean()), "paired_frac_mean": float(g.paired_frac.mean()),
            "n_stacks_mean": float(g.n_stacks.mean()), "stack_dG_mean": float(g.stack_dG_mean.mean()),
            "gc_stack_frac_mean": float(g.gc_stack_frac.mean()),
        })
    sdf = pd.DataFrame(summ).sort_values(["temp", "source"])
    sdf.to_csv(OUT / "mfe_summary.csv", index=False)

    # within-source association of MFE with complementary 6-mer pairs (37 C)
    assoc = {}
    for name, g in df[df.temp == 37.0].groupby("source"):
        x, y = g.cp_k6.to_numpy(float), g.mfe.to_numpy()
        assoc[name] = float(np.corrcoef(x, y)[0, 1])
    # between-source: does cp_k6 mean explain median MFE?
    kmer_rows = []
    for name in src:
        row = {"source": name}
        for k in KS:
            row[f"obs_k{k}"] = comp[name][f"obs_comp_pairs_k{k}"]
            if name in th:
                row[f"theory_k{k}"] = th[name][f"exp_comp_pairs_k{k}"]
        kmer_rows.append(row)
    pd.DataFrame(kmer_rows).to_csv(OUT / "comp_kmer_pairs.csv", index=False)
    pd.DataFrame({n: comp[n]["dinuc"] for n in src}).to_csv(OUT / "dinucleotide_freq.csv")

    out = {
        "settings": {"L": L, "N": N, "seed": SEED, "temps": TEMPS, "min_loop": MIN_LOOP,
                     "vienna_params": "DNA Mathews 2004, default salt (1.021 M Na+), dangles=2",
                     "nn_table": "SantaLucia 1998 unified dG37"},
        "theory": {k: {kk: vv for kk, vv in v.items() if kk != "dinuc"} for k, v in th.items()},
        "composition": {k: {kk: vv for kk, vv in v.items() if kk != "dinuc"} for k, v in comp.items()},
        "corr_mfe_vs_cp6_within_source_37C": assoc,
        "with_primers_37C": {k: {"median": float(np.median(v)), "ci": boot_ci(np.array(v))} for k, v in pres.items()},
    }
    (OUT / "summary.json").write_text(json.dumps(out, indent=2))

    # figures
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    for t, mk in zip(TEMPS, ("o", "s")):
        d = sdf[(sdf.temp == t) & sdf.p_repeat.notna()].sort_values("p_repeat")
        ax.errorbar(d.p_repeat, d.mfe_median, yerr=[d.mfe_median - d.mfe_median_ci_lo, d.mfe_median_ci_hi - d.mfe_median],
                    marker=mk, capsize=3, label=f"{t:.0f} °C (Markov sweep)")
        for name, col in (("iid_hp3", "C2"), ("shiftban", "C3")):
            r = sdf[(sdf.temp == t) & (sdf.source == name)].iloc[0]
            ax.errorbar([0.25 if name == "iid_hp3" else 0.0], [r.mfe_median],
                        yerr=[[r.mfe_median - r.mfe_median_ci_lo], [r.mfe_median_ci_hi - r.mfe_median]],
                        marker="D" if name == "iid_hp3" else "^", color=col, capsize=3, linestyle="none",
                        label=f"{name} {t:.0f} °C" if t == 37.0 else None)
    ax.set_xlabel("P(next base = previous base)")
    ax.set_ylabel("median MFE, 160 nt (kcal/mol), 95% bootstrap CI")
    ax.set_title("MFE vs repeat probability")
    ax.legend(fontsize=7)
    ax = axes[1]
    names = ["rot", "shiftban", "p0.05", "p0.10", "p0.15", "iid_hp3", "iid"]
    for k in KS:
        ax.plot(range(len(names)), [comp[n][f"obs_comp_pairs_k{k}"] for n in names], marker="o", label=f"k={k} observed")
    ax.set_xticks(range(len(names)), names, rotation=30)
    ax.set_yscale("log")
    ax.set_ylabel("complementary k-mer pairs / seq (loop ≥ 3)")
    ax.set_title("Hairpin-capable complementary k-mer pairs")
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"composition_mfe.{ext}", dpi=150)
    print(json.dumps(out, indent=1))
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
