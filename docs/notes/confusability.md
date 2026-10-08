# Nanopore confusability: definitions, provenance, findings (Phase 2)

Code: `dnastore/kmer.py`. Report: `scripts/confusability_report.py` → `results/phase2/confusability/`. Tests: `tests/test_kmer.py`.

## Provenance
`scripts/fetch_kmer_models.py` downloads from `nanoporetech/kmer_models` at pinned commit `4e56daed7fbb79b538f58e41262d5c54b07356ea` (master on 2026-10-01; last commit 2023-06-06). Provenance is recorded in `data/kmer_models/PROVENANCE.json`, and `load_model` re-verifies the SHA-256 on every load.

| name | file in repo | k | units | sha256 |
|---|---|---|---|---|
| `r10.4.1_9mer` | `dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt` | 9 | standardised (mean ≈ 0, sd ≈ 1) | `e5e15451…cf23c` |
| `r9.4.1_6mer` | `legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model` | 6 | pA (`level_mean`, `level_stdv`) | `5733049b…c16e` |

**Units caveat.** The R10.4.1 table has no pA scale and no per-k-mer std, so its δ is in standardised units and cannot be compared numerically with R9 δ in pA. ONT labels the legacy directory "r9.4". We call it r9.4.1 because it is the template model normally used for R9.4.1 flow cells.

**Orientation assumption.** Table k-mers are read 5′→3′ along the basecalled strand. Both strands of a dsDNA oligo are sequenced, so oligo metrics cover seq and revcomp(seq).

## Definitions
- **k-mer code:** A=0, C=1, G=2, T=3, MSB = first base. Rolling update: `code = ((code << 2) | b) & (4**k - 1)`.
- **RC lookup:** `rc_code_table(k)[c]` is the code of revcomp(k-mer c). Conf on the RC strand is `table[rc_code[codes]]`, so the beam never builds the RC string.
- **Confusability:** conf(x) = |{y : d_H(x,y) = 1, |ℓ(x) − ℓ(y)| < δ(x,y)}| / (3k), which lies in [0, 1].
  - R10.4.1: δ is constant, in std units (default 0.05, fixed a priori in DESIGN.md D1).
  - R9.4.1: δ(x,y) = c·√((σx² + σy²)/2), the pooled RMS of the two k-mers' `level_stdv` (default c = 1). This choice makes the relation symmetric, and it is the natural scale for separating two Gaussians of unequal width. Using σx alone would make "x confusable with y" asymmetric.
- **Oligo score** (`oligo_confusability`): mean and max over all k-mers of both strands, plus per-strand values.
- **Adjacent contrast** (`adjacent_contrast`, Whritenour, Civelek & Farnoud 2025, doi:10.1038/s41598-025-08531-z):
  - Definition: |μ(x_i) − μ(x_{i+1})| over consecutive overlapping k-mers.
  - What we matched: their constraint is a strict > τ, with μ from the ONT R9.4.1 6-mer means (pA) and τ ∈ {1…8} pA. We report min, mean, and `frac_below_tau` (the fraction of transitions with |Δμ| ≤ τ, i.e. violating their constraint), with τ = 4 pA by default for R9.
  - Our additions: (a) the RC strand, which they did not consider; (b) an R10 variant with τ = 0.5 std units, which is an **arbitrary placeholder** because there is no published τ for R10.

## Findings (1000 oligos per kind, seed 0)
Each oligo is 160 nt plus the configs/primers.yaml primers (200 nt total). Bodies are either i.i.d. quaternary (whitened) or rotating ternary chained from F[−1].

**Per-k-mer table distributions** (`kmer_table_stats.csv`):

| table | mean conf | median | p95 | frac 0 |
|---|---|---|---|---|
| R10 δ = 0.01 | 0.055 | 0.037 | 0.148 | 0.254 |
| R10 δ = 0.025 | 0.132 | 0.111 | 0.259 | 0.049 |
| **R10 δ = 0.05** | **0.244** | 0.222 | 0.444 | 0.005 |
| R10 δ = 0.1 | 0.423 | 0.407 | 0.630 | 0.000 |
| R10 δ = 0.2 | 0.646 | 0.667 | 0.815 | 0 |
| R9 c = 0.5 | 0.216 | 0.222 | 0.389 | 0.022 |
| R9 c = 1 | 0.404 | 0.389 | 0.611 | 0 |
| R9 c = 2 | 0.597 | 0.611 | 0.778 | 0 |

**Per oligo** (means over oligos):

| | R10 δ0.05 mean | R10 δ0.05 max | R9 c1 mean | R9 adjacent ≤ 4 pA fraction |
|---|---|---|---|---|
| whitened quaternary | 0.2419 | 0.579 | 0.405 | 0.158 |
| rotating ternary | 0.2425 | 0.572 | 0.433 | 0.104 |

1. **The per-oligo mean confusability barely varies.** Its std is 0.0045 on a mean of 0.242 (R10, δ = 0.05) and the full range across 1000 oligos is 0.226–0.257. Averaging over ~380 k-mers washes it out. The sparse steering bases can move the *mean* only a little. The beam should target the upper tail (max, or a sum of conf above a quantile) if the confusability objective is to have measurable effect. I'm recording this as a design input, not a decision.
2. **Alphabet does not matter for R10 confusability.** Means are 0.2419 vs 0.2425. Under R9, rotating ternary is slightly *worse* by our metric (0.433 vs 0.405) but **better** by Whritenour's: 10.4% vs 15.8% of transitions ≤ 4 pA, presumably because homopolymer steps produce small level changes.
3. **The two models disagree.** Per-oligo R10 vs R9 confusability has Spearman ρ = 0.10 (0.11–0.12 within each kind). "Confusability" is therefore model-specific, and claims must name the pore and model.
4. **Our metric and Whritenour's measure different things.**
   - Under R10, ρ is between −0.03 and 0.05 against both adjacent-contrast variants.
   - Under R9, conf correlates *positively* with the mean adjacent contrast: ρ = 0.90 over both kinds, 0.63–0.67 within each kind. Oligos with more substitution-confusable k-mers tend to have *larger* level steps. If both objectives were optimised together under R9 they could conflict. This is worth reporting.
5. **δ sensitivity is large.** The ranking of oligos at δ = 0.01 vs δ = 0.2 correlates at only ρ = 0.15 (0.05 vs 0.1: ρ = 0.67). This is why the δ sweep is mandatory: conclusions at one δ need not transfer.

## Caches
Tables are cached as `.npy` files in `data/kmer_models/cache/`, keyed by model, SHA-256 prefix, and parameter. The R10 9-mer table takes about 1–2 s to build (4^9 × 27 lookups, vectorised).
