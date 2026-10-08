# Structure metrics: definitions, validation, thresholds (Phase 2)

Code: `dnastore/structure.py`. Validation: `scripts/validate_structure_proxy.py` (results in
`results/phase2/structure/`: `proxy_validation.csv`, `proxy_validation_summary.json`,
`salt_sensitivity.txt`). Reproduce with
`.venv/bin/python scripts/validate_structure_proxy.py` (≈14 min on 8 cores), or
`--analyze-only` to reuse the CSV.

## 1. ViennaRNA conditions (2.7.2, verified)
- DNA parameters: `RNA.params_load_DNA_Mathews2004()`, loaded once per process.
- **Temperature:** `md.temperature` (°C). It changes energies as expected. For a test hairpin at 0.05 M: −21.19 kcal/mol at 37 °C, −12.36 kcal/mol at 60 °C.
- **Salt:** `md.salt` = **monovalent** cation concentration in mol/L (default 1.021 M). It changes energies: the same hairpin at 37 °C gives −24.80 kcal/mol at 1.021 M and −21.19 at 0.05 M. Related fields: `saltDPXInit`, `saltDPXInitFact`, `saltMLLower`, `saltMLUpper`, `helical_rise`, `backbone_length`.
- **No Mg²⁺ model.** The model-details object has no divalent-cation field. PCR buffers contain about 1.5–2.5 mM Mg²⁺, which stabilises structure strongly.
- **Defaults:** T = 60 °C (primary; PCR annealing) and 37 °C; salt = 0.05 M monovalent (the KCl in a standard PCR buffer, with Mg²⁺ ignored).
- **Salt sensitivity (important).** One way to fold Mg²⁺ in is a Na⁺-equivalent: Na_eq = Na⁺ + 120·√(Mg²⁺ − dNTP) in mM, attributed to von Ahsen et al. 2001 (*Clin. Chem.*). I quoted this formula from memory and **did not verify it**. It gives about 0.19 M for 50 mM K⁺ + 1.5 mM Mg²⁺ − 0.2 mM dNTP. On 101 oligos at 60 °C, median min-anchor accessibility is **0.710 at 0.05 M vs 0.462 at 0.19 M**, with Spearman 0.714 between the two. Absolute accessibility values and thresholds therefore depend strongly on the ion model. Rankings are only moderately stable. Report this sensitivity alongside any threshold-based result.

## 2. Primary metric: primer-site accessibility
Top strand T = F + body + rc(R), length L. Bottom strand B = rc(T) = R + rc(body) + rc(F).
- Primer R anneals to T at `[L−r, L)` (site = rc(R)). Primer F anneals to B at `[L−f, L)` (site = rc(F)).
- `P_u(i)` = 1 − Σ_j P(i·j), from the full McCaskill partition function (`fold_compound.pf()` + `bpp()`, with `exp_params_rescale(mfe)`).
- `site_access_X` = mean `P_u` over X's site. `min_site_access` = the minimum over F and R.
- **3′ anchor.** Primer and site are antiparallel, so the primer's 3′-terminal nucleotide pairs with the **5′-most** site position. Polymerase extends from there, towards the template's 5′ end.
  - `anchor_access_R` = mean `P_u` over `T[L−r : L−r+a]`.
  - `anchor_access_F` = mean `P_u` over `B[L−f : L−f+a]`, with a = 8.
  - `min_anchor_access` = the minimum over F and R. **This is the primary structure metric.**
- `structure_report` gives all of these plus MFE (forward and RC) at every configured temperature (keys suffixed `_T60`, `_T37`).

**Method choice.** For 200 nt, global pf + bpp takes ≈415 ms (one strand, uncontended). The windowed `probs_window` (RNAplfold-equivalent) at W = L = 200 is identical (max |Δ| = 4e-15) but slower (≈784 ms). With a 60-nt bp span it is faster but wrong (max |Δ| = 0.78). The global pf is used. Under 8-process load, a full accessibility evaluation (2 strands) took 2.37 s per oligo per temperature.

**Primer sanity check.** In isolation, rc(F) anchor P_u ≈ 1.0 and rc(R) anchor ≈ 0.84–0.99. In R + A₄₀/T₄₀ + rc(F) contexts, the anchors are ≥ 0.95. The primers are not intrinsically structured. **Random bodies** reduce the per-site medians to ≈0.82–0.86 (R anchor) and ≈0.78–0.82 (F anchor) at 60 °C.

## 3. Proxy validation (held-out half, n = 505 of 1010)
Oligo set (seed 20261001): 450 `quat` (15-nt rotating index + 145-nt i.i.d. quaternary body), 450 `tern` (160-nt rotating ternary), and 110 `adv` (quaternary body with an inserted 6–14-nt substring of R or rc(F)). Weights for combined proxies (λ) were chosen on the even-index half only. The ρ values below are on the odd half. 95% CIs come from 2000 bootstrap resamples. The target is `min_anchor_access_T60` unless noted. Negative ρ is expected for "more structure" scores.

| Proxy | Cost (ms/oligo) | ρ vs min_anchor_T60 [95% CI] | ρ vs min_site_T60 | ρ within quat / tern |
|---|---|---|---|---|
| **Tier 1, incremental k-mer** `combo_k4-6` (self + 8·site, k = 4,5,6) | 5.1 | **−0.373** [−0.450, −0.288] | −0.387 | −0.29 / −0.25 |
| `site_k4-6` alone | (part of above) | −0.304 [−0.380, −0.217] | −0.384 | −0.29 / −0.20 |
| `self_k6-7` alone | 3.7 | −0.206 [−0.294, −0.122] | −0.151 | −0.07 / −0.14 |
| **Tier 2, local pf, W = 60** | 48 | **+0.703** [+0.629, +0.771] | +0.612 | +0.86 / +0.69 |
| **Tier 2, local pf, W = 80** | 119 | **+0.778** [+0.710, +0.839] | +0.681 | +0.92 / +0.77 |
| Tier 2, local pf, W = 100 | 223 | +0.851 [+0.791, +0.900] | +0.754 | +0.95 / +0.84 |
| Limited-span MFE (span 30 / 50), reference only | 72 / 212 | +0.099 / +0.117 | ≈0 | — |
| Global MFE min(fwd, rc) at 60 °C | 796 | +0.207 [+0.120, +0.292] | +0.100 | +0.21 / +0.15 |
| Global MFE at 37 °C | 793 | +0.272 [+0.193, +0.349] | +0.171 | +0.37 / +0.20 |
| GC fraction | ~0 | −0.151 | −0.139 | — |

Findings:
1. **Global MFE is a poor surrogate for primer-site accessibility** (|ρ| ≈ 0.2–0.27). This supports making accessibility the primary structure metric.
2. **The cheap incremental k-mer proxy is weak** (|ρ| ≈ 0.37, below the 0.6 bar). Most of its signal comes from the primer-site complementarity term. It is only good enough as a soft tie-breaker.
3. **Local partition functions over the terminal W nt work well** (ρ = 0.70–0.85). The F-site term depends only on the **first** W nt of T, so a left-to-right beam can evaluate it exactly once position W is reached and prune there. The R-site term depends only on the last W nt and is evaluated at the end. W = 80 (≈60 ms per strand) is the recommended tier-2 default.
4. ρ is lower within `tern` than within `quat` for every proxy. Rotating-ternary bodies make long-range site sequestration more common, which the local window cannot see.

**Recommended tiering for the beam:** tier 1 (incremental k-mer proxy, soft cost) at every step; tier 2 (local pf W = 80; F site at position 80, R site at completion) on beam survivors; tier 3 (full `accessibility` at 60 °C, both strands) for final rescoring and reporting.

## 4. Chosen incremental proxy: exact definition
Parameters: `ks` (default (5, 6); the validation's best was (4, 5, 6) with `w_site` = 8), `min_loop` = 3, `w_self`, `w_site`.
`stem_weight(z)` = −Σ NN ΔG°37 over the k−1 Watson–Crick stacks of z (SantaLucia 1998 unified values), which is > 0.
When base i is appended, for each k with i ≥ k−1, let z = oligo[i−k+1 .. i]:
1. First, if e = i − k − min_loop ≥ k−1, add the k-mer ending at e (oligo[e−k+1 .. e]) to the count table `C_k`.
2. Self term: `w_self · stem_weight(z) · C_k[rc(z)]`. This counts earlier complementary k-mers separated by a loop of at least `min_loop`.
3. Site term: `w_site · stem_weight(z) · [z ∈ S_k]`, where `S_k` = k-mers of R ∪ k-mers of rc(F), precomputed from the primers. A T k-mer pairs with the R site iff it is a k-mer of R. It pairs with the F site on B iff it is a k-mer of rc(F). The term is strand-symmetric.

The increment at i is the sum over k of (2) + (3). `proxy_score` = Σ increments (`proxy_increments` returns the list; pass `primers=(F, R)` when scoring a prefix). Each step is O(|ks|). A test checks the prefix property: appending bases never changes earlier increments.

## 5. Threshold calibration (from accessibility, not MFE percentiles)
**Model (pessimistic, explicit):** per-cycle amplification efficiency is proportional to equilibrium 3′-anchor accessibility a (E = a·E_max, E_max = 1). After n cycles, an oligo's yield relative to a reference with accessibility a_ref is ((1+a)/(1+a_ref))ⁿ. To be at most `fold`× under-represented, a ≥ (1 + a_ref)·fold^(−1/n) − 1.

- **Absolute version (a_ref = 1, a fully accessible site):** n = 25, fold = 2 gives a* = 0.945. **It rejects 100% of oligos, including all unconstrained random ones** (the median min-anchor accessibility is 0.72). Either the pessimistic model overstates the effect (primer excess may kinetically outcompete partial intramolecular structure; equilibrium P_u at 0.05 M with Mg²⁺ ignored is itself uncertain), or every random oligo is heavily disadvantaged against an ideal one. Both are untested, so the absolute version is not used as a hard constraint.
- **Proposed default (relative version):** a_ref = median min-anchor accessibility of the unconstrained bodies (quat + tern) = 0.717. With n = 25 and fold = 2: **a* = 0.686 on `min_anchor_access_T60`**. In words: reject an oligo predicted, under the pessimistic model, to be ≥ 2× under-represented after 25 cycles relative to a typical unconstrained oligo. It was fixed before any codec results existed.

Fraction of oligos failing (min_anchor_access_T60 < threshold):

| Rule | a* | quat | tern | adv |
|---|---|---|---|---|
| rel n = 25, fold 2 (**proposed**) | 0.686 | 0.371 | 0.420 | 0.736 |
| rel n = 25, fold 4 | 0.636 | 0.309 | 0.322 | 0.691 |
| rel n = 20, fold 2 | — | 0.356 | 0.389 | 0.727 |
| rel n = 30, fold 2 | — | 0.387 | 0.449 | 0.755 |
| abs n = 25, fold 2 | 0.945 | 1.000 | 1.000 | 1.000 |
| abs n = 25, fold 4 | — | 0.887 | 0.931 | 0.973 |

Distributions of min_anchor_access_T60 (p1 / p5 / p25 / median): quat 0.24 / 0.37 / 0.60 / 0.75; tern 0.28 / 0.41 / 0.60 / 0.71; adv 0.02 / 0.04 / 0.22 / 0.41.

**Recommendation:** treat a* as a **sweep parameter** (rows above) rather than a fixed hard constraint until the coverage model (Phase 2 item 4, from published PCR-bias data) links accessibility to measured copy-number bias. At the proposed a* ≈ 0.69, about 40% of unconstrained oligos would fail. That is a large cost for the steering beam or seed rerolls and must be reported as such.

## 6. Limitations
- No Mg²⁺ in ViennaRNA. Results are sensitive to the ion model (§1).
- Accessibility is an equilibrium, single-molecule ssDNA property. In PCR, primer binding competes kinetically with intramolecular folding at a high primer concentration, the template is single-stranded only transiently after denaturation, and the polymerase can open weak structure. Equilibrium P_u is a **proxy**.
- Mean P_u over the anchor positions is not the joint probability that all of them are unpaired. The joint probability would be stricter; Vienna's `probs_window` with `ulength` = 8 could compute it at extra cost.
- The efficiency-vs-accessibility link (§5) is an assumed model, not a fitted one.
- DNA Mathews 2004 parameters are less well validated than RNA Turner parameters, particularly for G·T and dangles.
- The validation set is synthetic (random bodies plus planted complementarity), not steering-codec outputs. Re-check ρ on actual steering-B outputs once they exist, because the beam may exploit proxy blind spots (Goodhart).
