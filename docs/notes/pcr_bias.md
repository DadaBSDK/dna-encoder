# PCR / synthesis coverage-bias model: sources, parameters, limitations

Module: `dnastore/coverage.py`. Tests: `tests/test_coverage.py`. Written 2026-10-01.

**Purpose.** squigulator + Dorado simulate *sequencing* only. A constraint-violating arm (whiten+RS) could pay a hidden cost in amplification yield. This model lets per-oligo copy number depend on GC content and primer-site accessibility. **It is run only as a sensitivity sweep over the presets below, never as a single "true" channel.**

## 1. Sources (every number used, with its location)

All numbers were read from the source text (PMC full text or the PDF), not from memory. Quotes are verbatim or near-verbatim from the fetched text.

| # | Source | Quantity | Value | Location | Conditions |
|---|---|---|---|---|---|
| S1 | Gimpel et al. 2023, "A digital twin for DNA data storage…", *Nat Commun*. [doi:10.1038/s41467-023-41729-1](https://doi.org/10.1038/s41467-023-41729-1) ([PMC10533828](https://pmc.ncbi.nlm.nih.gov/articles/PMC10533828/)) | Synthesis coverage, lognormal shape σ | Twist (material deposition): **0.27** unconstrained, **0.30** GC-constrained. CustomArray/GenScript (electrochemical): **1.30** unconstrained, **0.58** GC-constrained | Fig. 2b and Results ("All pools fit a lognormal distribution") | Pools from the authors' experiments |
| S1 | same | Relative PCR efficiency (1+ε_i)/(1+ε̄), sd | **0.0051** (unconstrained), **0.0048** (GC-constrained), normally distributed. Literature datasets: **0.0058–0.012** | Fig. 3d and caption | Serial PCR: 15 cycles per round, six rounds |
| S1 | same | GC dependence of PCR bias | "the broadness of the efficiency distribution does not appear to directly depend on GC constraints" | Fig. 3d text | — |
| S1 | same | Stochastic PCR | "branching binomial process, based on oligonucleotide count and the sequence's amplification efficiency" | Methods | — |
| S2 | Chen et al. 2020, "Quantifying molecular bias in DNA data storage", *Nat Commun* 11:3264. [doi:10.1038/s41467-020-16958-3](https://doi.org/10.1038/s41467-020-16958-3) ([PMC7324401](https://pmc.ncbi.nlm.nih.gov/articles/PMC7324401/)) | PCR model | n_{j+1} = n_j + B(n_j, P), "P is the probability of a successful amplification"; **P = 0.95** in the Fig. 5b model. **The paper does not say whether P was fitted or assumed.** | Methods; Fig. 5b | — |
| S2 | same | GC effect | "practically unimportant (the slope of the linear fit was <0.01)", attributed to the "relatively short (150-nt)" oligos and KAPA HiFi | Fig. 4c,d text | 150-nt oligos, GC 25–75%, KAPA HiFi |
| S2 | same | Coverage c.v. | 0.41 (before) and 0.45 (after PCR) | Fig. 4a | Optimised synthesis process |
| S2 | same | Dominant PCR bias | Stochasticity at low copy number ("the lower oligo copy numbers were, the greater the PCR stochastic bias was") | Fig. 5d caption | — |
| S3 | Gimpel et al. 2025, "Predicting sequence-specific amplification efficiency in multi-template PCR with deep learning", *Nat Commun*. [doi:10.1038/s41467-025-64221-4](https://doi.org/10.1038/s41467-025-64221-4) ([PMC12533003](https://pmc.ncbi.nlm.nih.gov/articles/PMC12533003/)) | Poor-amplifier tail | "around **2%** of the pool" with "efficiencies as low as **80%** relative to the population mean (equivalent to a halving in relative abundance every 3 cycles)" | Fig. 2c | KAPA SYBR FAST; 149-nt oligos; 6 × 15 cycles |
| S3 | same | GC as predictor | "cannot be explained by the base composition or GC content alone" (regression close to a random classifier) | Results, "Positional sequence information is critical" | — |
| S3 | same | Mechanism | "adapter-mediated self-priming": motifs (CGTG) complementary to the 5′ end of the 5′ adapter form hairpins next to the primer site; predicted hairpin Tm 45–53 °C, "competitive to primer annealing" | Fig. 4c,f | — |
| S3 | same | Motif effect size | Inserting TCGTGT "led to a mean decrease in PCR efficiency of **4.8±2.4%**" | Fig. 6c | External-lab validation |
| S4 | Heckel, Mikutis & Grass 2019, *Sci Rep* 9:9663. [doi:10.1038/s41598-019-45832-6](https://doi.org/10.1038/s41598-019-45832-6) (read from [arXiv:1803.03322v1](https://arxiv.org/abs/1803.03322), Sec. 4) | PCR efficiency | Factor E ~ Gaussian(**1.85**, **0.07**), 22 and 60 cycles. **This is a simulation assumption**, citing [Rui+09; Pan+14; War+97; Cal+07], not a measurement | Sec. 4, Fig. 8b,c | Illustrative |
| S5 | Aird et al. 2011, *Genome Biol* 12:R18. [doi:10.1186/gb-2011-12-2-r18](https://doi.org/10.1186/gb-2011-12-2-r18) | GC bias | "as few as ten PCR cycles … (Phusion HF …) … depleted loci with a GC content > **65%** to about a **hundredth** of the mid-GC reference loci … Amplicons < **12%** GC were diminished to approximately **one-tenth** … plateau … ranged from **11% to 56%** GC" (plateau = relative abundance ≥ 0.7) | Results, Fig. 1e | **Genomic** fragment libraries, standard Illumina Phusion protocol. The optimised protocol (betaine, longer denaturation) largely removed the high-GC bias (Fig. 5) |

## 2. Model

For oligo *i* with GC fraction *g_i* and primer-site accessibility *a_i* ∈ [0, 1] (computed by the caller; e.g. the minimum unpaired probability over primer sites):

```
s_i  = n0 · LogNormal(−σ_syn²/2, σ_syn)                      # synthesis, mean n0            (S1)
r_i  = (1 + N(0, sd_rel)) · g_gc(g_i) · g_acc(a_i) · tail_i  # relative per-cycle efficiency (S1 definition)
e_i  = clip(r_i · (1 + ē) − 1, 0, 1)
deterministic:  N_i = s_i · (1 + e_i)^c
stochastic:     n ← Poisson(s_i); repeat c times: n ← n + Binomial(n, e_i)   (S2; expectation used once n ≥ 1e4)
g_gc  = 1                                   (DNA-storage evidence: S1, S2, S3)
      | Aird curve (out-of-domain bound):   per-cycle factor 0.1^(1/10) for GC ≤ 0.11, 1 on [0.12, 0.56],
                                            0.01^(1/10) for GC ≥ 0.65, log-linear between (ASSUMED shape between anchors)
g_acc = 1 − κ · clip((a0 − a_i)/a0, 0, 1)   (ASSUMED form; κ magnitudes anchored to S3)
tail_i = tail_rel_eff with prob. tail_frac, else 1   (S3)
```

## 3. Parameter table and presets (`SENSITIVITY_GRID`)

| Parameter | Value(s) | Source | Status |
|---|---|---|---|
| σ_syn | 0.27 (Twist), 1.30 (CustomArray) | S1 Fig. 2b | literature-measured |
| sd_rel | 0.0051 (measured), 0.012 (upper end across datasets) | S1 Fig. 3d | literature-measured |
| ē | 0.95 | S2 Fig. 5b (origin of P not stated) | literature value, provenance weak |
| ē, sd_rel (Heckel) | 0.85, 0.07/1.85 = 0.038, 22 cycles | S4 Sec. 4 | modelling assumption in source |
| c (cycles) | 15 (one round) | S1 Methods | literature protocol |
| GC curve | none (default) | S1, S2, S3 | literature: no practical GC effect for ~150-nt storage oligos with KAPA |
| GC curve | Aird anchors | S5 | literature, but **out of domain** (genomic libraries, Phusion) → pessimistic bound only |
| tail | 2% at relative efficiency 0.8 | S3 Fig. 2c ("as low as 80%": the lower bound is used, i.e. pessimistic) | literature-measured, applied at random (not sequence-linked) |
| κ (accessibility penalty) | 0.048, 0.20 | magnitudes from S3 (Fig. 6c motif effect; Fig. 2c tail) | **ASSUMED: no published accessibility→efficiency function** |
| a0 (accessibility threshold) | 0.5 | none | **ASSUMED** |
| n0 | 100 molecules per oligo | none (sampling choice) | free parameter; affects only stochastic dropout. Sweep it |

Presets: `uniform`, `twist_synth_only`, `twist_pcr_measured`, `pcr_spread_lit_upper`, `customarray_pcr`, `poor_tail`, `heckel_illustrative`, `aird_gc_pessimistic`, `acc_assumed_mild`, `acc_assumed_strong`. Each preset's `sources` dict tags its non-default values. Presets with κ > 0 carry the "ASSUMED" tag, and a test enforces this.

Smoke check (20k oligos, GC ~ N(0.5, 0.035), accessibility ~ U(0, 1); weights normalised to mean 1): sd(log w) is 0.27 (synth only), 0.30 (measured PCR), 0.34 (upper spread), 0.58 (poor tail; 2.2% of oligos below 0.1×), 0.42 (Aird), 0.38 / 1.10 (accessibility mild / strong). This is a model output, not a result.

## 4. Limitations (state these in the paper)

1. **squigulator + Dorado model sequencing only.** This module models amplification *yield* only: copy number per oligo. It does **not** model synthesis errors, PCR substitutions (S1 reports ~1.09e-4 per nt per cycle; not used here), chimeras, mispriming products, length truncation, or storage decay.
2. **No published quantitative accessibility → efficiency function was found.** S3 shows the *direction*: hairpins between the template near the primer site and the adapter lower efficiency. It gives effect sizes for specific motifs, but no function of a computed accessibility or unpaired probability. The dependency is therefore an **assumed knob** whose magnitude is merely bounded by S3. Any accessibility-driven result is conditional on that assumption and must be reported as such.
3. **The GC effect is not supported for storage oligos under KAPA chemistry** (S1, S2, S3). The Aird curve comes from genomic fragment libraries with Phusion under a non-optimised protocol, and uses amplicon GC rather than whole-oligo GC. It is a pessimistic out-of-domain bound. Also, whitened 200-nt oligos almost never leave the 40–60% global GC range, which lies inside Aird's plateau (except > 56%). So the GC term mostly matters for oligos with high local GC, which this model, being based on whole-oligo GC, does not resolve.
4. **The poor-amplifier tail is applied at random** with the measured prevalence. In reality it is sequence-linked (S3's CNN). A sequence-linked version would need S3's model or data, which we have not integrated.
5. **P = 0.95 (S2) has no stated provenance.** The deterministic-limit results scale with (1+ē)^c only through relative efficiencies. ē matters mainly for stochastic low-copy dropout. Sweep it if conclusions depend on it.
6. Lognormal synthesis spread is applied independently of sequence. S1/S2 report spatial (chip-location) effects, which are not modelled.
