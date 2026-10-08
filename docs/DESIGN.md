# DESIGN: Structure- and Readout-Aware DNA Storage Codec

Phase 0 deliverable (2026-10-01), revised after review (Revision 1, same date). **Part 0 supersedes anything in Parts 1–2 that conflicts with it.** Parts 1–2 are kept as the record of the original critique.

See `RELATED_WORK.md` first. Its overlap warning is the reason for the framing in Part 0.

Implementation update (2026-10-08): the matched-budget runner described in §2.6 is now
available as `dnastore benchmark`. See [BENCHMARK.md](BENCHMARK.md) for exact budget
matching and statistical interpretation, and [PROJECT_STATUS.md](PROJECT_STATUS.md)
for verified runs and remaining research work. The complete publication study is still
pending; smoke and regression experiments do not resolve the research question.

---

## Part 0 — Research question and evaluation plan (Revision 1, binding)

### 0.1 Core question
> **Do structure- and readout-aware sequence constraints pay for the nucleotides they cost, compared with spending the same nucleotides on error correction?**

The steering-base codec is the **instrument** for answering this. The steering period P sets how many nucleotides go to constraints: 1/(P+1) of body positions. P = ∞ spends none.

### 0.2 Comparison arms (all at equal total synthesised nt per file)
| Arm | Role |
|---|---|
| **steering-B** (rotating ternary + steering), P sweep | **Primary instrument.** Homopolymers are handled by construction (run ≤ 2), so every steering base is free for GC, structure, and confusability. |
| **whiten+RS** (whitened, unconstrained payload; nt saved from steering go to extra outer-RS parity) | **Primary comparison arm** (the Weindel et al. "embrace errors" hypothesis). Matched to each steering-B configuration by total nt. |
| steering-A (quaternary + whitening + seeds), P sweep | Secondary arm. Steering helps homopolymers only for P ≤ 3 (F4); seeds do that work otherwise. |
| goldman (rotating ternary, P = ∞, shared ECC) | Baseline. Equal to steering-B at P = ∞. |
| naive2bit | Baseline: no constraints. |
| fountain (own LT + screening) | Baseline: screening-by-seed with rateless redundancy. |

"Equal total nt" means: fix the oligo length and pool size (oligo count) for a given file. Each arm divides its body nucleotides among payload, constraints, and parity differently. A constraint is judged worth its cost only if the steering arm beats the matched whiten+RS arm on the primary metrics (0.4).

### 0.3 Evaluation channels
1. **Main read-level evaluation: identity-matched Badread.** Badread v0.4.1's nanopore2023 model is calibrated to public R10.4.1 storage reads (mean identity 98.48%, SD 1.17%). The source data are plasmid fragments, not a 200-nt oligo pool, and Badread has no read-to-read systematic errors, so this remains an approximation.
2. **Controlled channel: uniform IDS.** Substitution, insertion, deletion, dropout, lognormal coverage, and 50% reverse-complement reads. Used for file recovery and minimum-depth sweeps.
3. **Stress channel: squigulator → Dorado.** Retained as a raw-signal/basecaller integration check. Calibration measured 9.12% per-base error versus 1.51% in real hac reads, plus a much larger systematic homopolymer floor; results must be labelled miscalibrated and cannot support the main readout claim. See `docs/notes/channel_calibration.md`.
4. **Confusability-weighted channel: sanity check only, labelled CIRCULAR** in every table and figure caption. It shows that the optimiser does what it is told. It is never evidence of better readout.

### 0.4 Metrics
- **Primary (continuous):** (a) post-consensus per-oligo error rate (edit distance of the consensus body vs. the true body ÷ body length, plus the fraction of oligos failing CRC), as a function of depth and error rate; (b) **minimum sequencing depth for 100% file recovery**, found by bisection over depth, with ≥ 20 channel seeds per probe (100 where affordable) and Wilson CIs on the recovery fraction at each probe.
- **Secondary:** binary file-recovery rate (Wilson 95% CI).
- **Descriptive:** densities (§2.5), constraint distributions, MFE, confusability (mean/max, δ sweep), **Whritenour et al. adjacent-level contrast on every codec**, runtimes, mass.

### 0.5 Structure proxy inside the beam (validated)
Full MFE (~85 ms per 200-nt fold) is used only for final rescoring. Inside the beam, a cheap proxy is used: candidates are (i) ViennaRNA folding with a limited base-pair span (`md.max_bp_span`, e.g. 30–50) on the local window, and (ii) a self-complementary stem scan (reverse-complement k-mer pairs, k = 4–6, weighted by stem length and GC). **Validation:** on ≥ 1000 oligos (whitened random bodies with our primers, plus steering outputs), report the Spearman ρ of each proxy against full MFE (fwd and RC) along with timing. The proxy actually used, its ρ, and its speed-up are reported in RESULTS.md. If ρ < 0.6, say so and reconsider.

### 0.6 Confusability
R10.4.1 9-mer (std. units) is primary, with **δ sweep** {0.01, 0.025, 0.05, 0.1, 0.2} (the 0.05 default is fixed a priori). R9.4.1 6-mer (pA, δ = c·level_stdv) is secondary. **Whritenour-style adjacent-level contrast** (min and mean over i of |ℓ(x_i) − ℓ(x_{i+1})|) is computed on all codecs.

### 0.7 Accepted resolutions
D1–D5 use the Phase 0 defaults (table at the end), plus the δ sweep. Header: R distinct, codec-independent header oligos (F6). The orientation marker is **dropped** (F8).

### 0.8 Design changes made in Phase 1 (flagged for review)
- **Shared, codec-independent index field.** Every oligo body starts with a fixed-width index field (rotating ternary, 15 trits ≈ 14.3 M indices; the first base is chained from the forward primer's last base). This deviates from the spec's "index uses the same constrained encoding". Reasons: (i) the decoder must group reads and find header oligos *before* it knows the codec; (ii) all codecs then pay an identical index cost (F11 fairness). The index is homopolymer-free by construction. Its GC and structure enter the steering beam as a fixed prefix.
- **Header oligos are not steered.** The header uses rotating ternary + whitening keyed by copy number (so the R copies differ completely), with no steering. It is a handful of oligos; their constraint values are screened and reported like any other oligo. Indices 0–63 are reserved for header fragments (index = fragment·R + copy).
- **Fixed oligo length across codecs.** Codecs whose body does not exactly fill the length (e.g. 2-bit mapping) pad with a fixed filler, and the padding is counted in density.

- **Header integrity.** There is no separate header CRC32. Each header fragment is protected by its oligo's CRC-16 (bound to the index), the parse checks magic and version, and the file-level CRC32 catches anything left. This saves a few bytes of header and changes §2.4 accordingly.
- **Custom GF(256) RS instead of reedsolo.** reedsolo is pure Python here (~6 ms encode / ~18 ms decode per codeword). Oligo dropouts erase the same row in every column of a block, so `ecc.py` implements RS as a Vandermonde parity-check code. It needs one erasure solve per block, vectorised over columns, plus single-bad-row location using spare syndromes. Measured on 100 KB: 0.09 s encode, 0.18 s decode at 5% dropout.

### 0.9 Phase 1 observations that affect later phases
- **Rotating ternary folds more stably than quaternary.** Random 160-nt sequences without primers have a median MFE of −26.6 kcal/mol for rotating ternary vs −17.4 for i.i.d. quaternary (150 each, ViennaRNA DNA Mathews 2004, 37 °C). Full 200-nt oligos with primers on 20 KB of random data, using min(fwd, RC): goldman has a median of −32.8 vs −22.6 for whitened 2-bit. So alphabet B, the primary arm, **starts from a worse structure baseline**, and its steering bases have more to fix. This should be reported, not hidden.
- **Baseline hard-constraint violation rates** (20 KB random bytes, hard = homopolymer, gc_global, gc_window(20 nt, 25–75%), primer ≤ 3 mm): naive2bit 85.1%; whitened 2-bit 88.5% (consistent with F4's ≈ 17% homopolymer pass rate); goldman 3.4% (all GC-window).

---

## Part 0b — Revision 2 (Phase 2 review), binding

### 0b.1 Candidate result: "homopolymer constraints carry a hidden structure cost"
> **Homopolymer constraints can carry a hidden structure cost, and the size depends on how the constraint is enforced.** Rotating ternary (run length 1) folds about 10 kcal/mol more stably than whitened quaternary: median MFE −27.1 vs −17.2 for 160 nt at 37 °C, −14.7 vs −8.8 at 60 °C, and −31.8 vs −22.9 for full 200-nt oligos with primers. The cause is lower entropy in the first-order chain. Banning one of four successors raises the transition collision probability from 1/4 to 1/3, which multiplies the expected number of hairpin-capable reverse-complementary k-mer pairs by (4/3)^(k−1): 4.2× for 6-mers, and theory matches observation. GC stacking is not the cause; `rot` has fewer GC stacks. A control that bans a different successor reproduces ~79% of the gap. Rejection-filtering to the conventional max run ≤ 3 has no measurable cost. **Implication:** alphabet B starts from a structural deficit that its steering bases must pay back. Any comparison of B against whitened quaternary has to charge B for it.
Source: `docs/notes/composition.md`, `results/phase2/composition/` (n = 2000 per source; a repeat-probability dose-response gives median MFE −17.2, −18.5, −20.7, −23.3, −27.1 for P(repeat) = 0.25, 0.15, 0.10, 0.05, 0). The claim is precise only for **run-length-1** enforcement. Refined wording: *homopolymer constraints enforced at every position carry a hidden structure cost; max-run ≤ 3 filtering does not.* The beam's stem proxy (complementary k-mer pairs) targets exactly this mechanism.

### 0b.2 Structure metrics
Explicit temperatures: report at **60 °C (PCR annealing, primary)** and **37 °C**, with ViennaRNA salt correction if supported. **Primary structure metric = primer-site and 3′-end accessibility** (unpaired probability from the partition function). Global MFE is secondary. Thresholds are calibrated from accessibility, not from MFE percentiles. The beam proxy is validated against the primary metric (Spearman, about 1000 oligos). Details are in `docs/notes/structure.md`; see §0b.7 for the results.

### 0b.3 3′-anchored primer matches are hard for ALL arms
Every oligo has a shared, codec-independent **seed field** (3 rotating-ternary nt, 27 seeds) after the index. The layout is `[F][index 15][seed 3][payload][rc(R)]`. Every arm (naive2bit, whiten+RS, goldman, steering-A/B, header oligos) rerolls its seed until there are no exact 3′-anchored primer 8-mer matches (F[−8:], R[−8:], rc(F)[:8], rc(R)[:8], true sites excluded). For baselines, seed 0 is the plain code and seeds > 0 apply keystream whitening. If all 27 seeds fail, the least-violating candidate is emitted and the violation is recorded. Header oligos, which every arm shares, reroll against *all* hard constraints. The whiten+RS arm otherwise keeps its violations: that is the hypothesis under test. The 3-nt seed field costs every arm the same.

### 0b.4 PCR/synthesis coverage model (optional, sensitivity only)
`dnastore/coverage.py`; sources and parameter table in `docs/notes/pcr_bias.md`. Literature-derived terms:
- synthesis lognormal σ (Gimpel 2023, Fig. 2b: 0.27 Twist, 1.30 CustomArray);
- PCR relative-efficiency spread sd_rel ≈ 0.005 (Gimpel 2023, Fig. 3d);
- binomial branching (Chen 2020);
- **no GC effect within storage-oligo ranges** (Chen 2020, Fig. 4; Gimpel 2023/2025);
- a 2% poor-amplifier tail at relative efficiency 0.8 (Gimpel 2025, Fig. 2c).

**Accessibility → efficiency has no published quantitative source.** It is an explicit knob labelled "ASSUMED" (κ ∈ {0.048, 0.20}; the magnitudes are borrowed from Gimpel 2025's motif and tail effects; threshold a0 = 0.5 is assumed). Aird 2011's GC curve is included only as an out-of-domain pessimistic bound.
**Limitation (stated in every results section that uses it):** squigulator + Dorado models sequencing only. This model adds amplification *yield* only. Neither models synthesis errors, PCR substitutions, chimeras, or decay. The whiten+RS arm's unconstrained oligos are therefore penalised only to the extent that these yield models allow, and with the literature finding no GC effect, that will likely be little. If whiten+RS wins, the conclusion must be qualified as "under sequencing-error and yield models; synthesis-side costs of unconstrained sequences are not modelled."

### 0b.5 RS correctness against a reference
The `ecc.py` convention was changed to match reedsolo (fcr = 0, prim 0x11D, α = 2, message first). `tests/test_ecc.py::test_rs_codewords_identical_to_reedsolo` checks that codewords are byte-identical to `reedsolo.RSCodec` on 80 random (k, p, message) draws.

### 0b.6 Note for Phase 3: index/seed substitutions
In rotating ternary, symbol *i* is the base *relative to base i−1*. **One substitution at position i corrupts two decoded symbols** (trits i and i+1), or is detected outright when it creates a repeat, which makes the step invalid. This applies to the 15-nt index and the 3-nt seed field. The grouping-error analysis in Phase 3 must account for it: misgrouping probability per read, mis-assignment to *another valid index* (caught by the index-bound CRC-16), and the edit-distance clustering fallback.

### 0b.7 Phase 2 results so far
- **Confusability** (`docs/notes/confusability.md`): tables pinned at commit `4e56dae`, sha256 recorded in `data/kmer_models/PROVENANCE.json`. Findings:
  - The per-oligo *mean* confusability barely varies (std 0.0045 on a mean of 0.242, R10 at δ = 0.05), so the beam penalises only the **upper tail** (excess above the 75th percentile of the per-k-mer table).
  - R10 and R9 confusability correlate at only ρ = 0.10, so any claim must name the pore model.
  - Our Hamming-neighbour metric and Whritenour's adjacent-contrast metric are nearly uncorrelated under R10 but correlated under R9 (ρ = 0.63–0.67 within an alphabet).
  - Rankings are δ-sensitive (ρ = 0.15 between δ = 0.01 and δ = 0.2), so the δ sweep is necessary.
- **Structure** (`docs/notes/structure.md`, `results/phase2/structure/`):
  - ViennaRNA 2.7.2 supports `md.temperature` and monovalent `md.salt` (mol/L; default 1.021) but **has no Mg²⁺ model**. Defaults are 60 °C (primary) and 37 °C, with 0.05 M monovalent salt. Salt matters: at 0.19 M (a rough Mg²⁺-equivalent; the conversion is unverified), median minimum anchor accessibility falls from 0.71 to 0.46, and the rank agreement with 0.05 M is only ρ = 0.71. **Absolute accessibility values depend on the ion model.**
  - The primary metric is accessibility of the primer sites and 3′ anchors on both strands (full partition function). The R site is the last r nt of the top strand. The F site is the last f nt of the bottom strand. The anchor is the first 8 positions of each site.
  - **Proxy validation** (1010 synthetic oligos, held-out half, Spearman vs min anchor accessibility at 60 °C):

    | Proxy | ρ |
    |---|---|
    | Tier-1 incremental proxy: NN-weighted self-complementary k-mer pairs + primer-site complementarity, k = (4, 5, 6), w_site = 8 | **−0.37** [−0.45, −0.29], **below the 0.6 bar** |
    | Tier-2 windowed partition function, W = 80 | **0.78** [0.71, 0.84] |
    | Tier-2 windowed partition function, W = 100 | 0.85 |
    | Global MFE | ≈ 0.2 |

    Global MFE is a poor surrogate for accessibility, which supports making accessibility primary. Proxies do worse inside rotating-ternary bodies.
  - **How the beam uses them:** the tier-1 proxy is a soft cost at every step. Its incremental implementation is tested to equal `structure.proxy_increments` exactly. The final `rescore_top` = 4 candidates are re-ranked by tier 2 (W = 80, 60 °C). Tier 3 (full accessibility at 60 and 37 °C) is used for reporting.
  - **Threshold:** an absolute, physically motivated threshold (≤ 2× under-representation over 25 cycles, assuming efficiency ∝ anchor accessibility) gives a* = 0.945. That rejects **100%** of oligos, random ones included (median 0.72). The relative version, referenced to the median of unconstrained bodies, gives **a* = 0.686**, which fails 37% of whitened-quaternary and 42% of rotating-ternary oligos. Both rest on an *assumed* efficiency ∝ accessibility link, which has no published quantitative support (§0b.4). **Decision (proposed): accessibility is a soft, rescoring objective, not a hard constraint.** The fraction of oligos with min anchor accessibility ≥ a* is reported for a* ∈ {0.6, 0.686, 0.8} as a sweep. Making it hard would require a reliable cheap check (tier 2 has ρ = 0.78) and a data-linked threshold; neither exists.
- **Phase 2 characterisation sweep** (`scripts/steering_sweep.py`; `results/phase2/sweep/`; 4 KB random bytes, seed 0, 200-nt oligos, data oligos only; one payload, so **descriptive, not a benchmark**):

  | arm | code b/nt | violating | mean tries | median min anchor access 60 °C | ≥ 0.686 | median MFE 37 °C fwd | stem proxy | conf max (R10) |
  |---|---|---|---|---|---|---|---|---|
  | naive2bit | 1.847 | 89% | 1.0 | 0.821 | 73% | −7.7 | 1451 | 0.580 |
  | whiten 2-bit (whiten+RS arm) | 1.847 | 91% | 1.0 | 0.798 | 70% | −7.2 | 1404 | 0.572 |
  | goldman (= B, P = ∞) | 1.463 | 3.4% | 1.0 | 0.769 | 68% | −14.7 | 1626 | 0.573 |
  | steering-B P = 4 | 1.127 | 0 | 1.0 | 0.828 | 83% | −7.1 | 359 | 0.569 |
  | steering-B P = 6 | 1.238 | 0 | 1.0 | 0.847 | 83% | −6.9 | 413 | 0.576 |
  | steering-B P = 8 | 1.295 | 0 | 1.0 | 0.840 | 80% | −7.9 | 522 | 0.570 |
  | steering-B P = 12 | 1.349 | 0 | 1.0 | 0.824 | 80% | −9.6 | 736 | 0.571 |
  | steering-A P = 4 | 1.463 | 0 | 1.7 | 0.865 | 88% | −3.9 | 479 | 0.577 |
  | steering-A P = 8 | 1.633 | 0 | 2.8 | 0.850 | 79% | −4.8 | 858 | 0.573 |
  | steering-A P = ∞ | 1.847 | 0 | 6.8 | 0.847 | 81% | −6.7 | 1348 | 0.570 |

  Readings (single payload; CIs come in Phase 4):
  1. **Steering-B recovers B's structural deficit.** Its median MFE at 37 °C goes from −15.0 at P = ∞ to −6.9 to −7.9 at P = 4–8, which is equal to or better than whitened quaternary. Anchor accessibility goes from 0.75 to 0.83–0.85. The price is density: 1.24–1.30 vs 1.85 b/nt for whiten+RS (−30%). Phase 4 asks whether that pays.
  2. **Steering-A is denser and structurally better than B** at equal or greater density: A at P = 4 matches goldman's density but has better accessibility (0.865 vs 0.769). It needs seed rerolls (mean 1.7–6.8, max 26), and the 27-seed field was never exhausted in this sample. *This contradicts the "B primary" premise on these metrics.* Flag for review.
  3. **A at P = ∞ is not an unconstrained baseline:** seed selection plus tier-2 rescoring is itself structure optimisation. The unconstrained arm is `naive2bit(whiten=True)`.
  4. **Byte quantisation:** P = 12 and P = 16 give identical density in both alphabets, because frames are whole bytes. P = 16 wastes steering freedom. Choose P values on the density lattice.
  5. **Confusability is NOT steered at the current weights.** conf max is ≈ 0.57 in every arm. A follow-up (800 B, B, P = 6) shows it *is* steerable when it is the only objective: tail excess drops from 4.42 to 1.80 (−59%), conf max from 0.576 to 0.498, and conf mean from 0.241 to 0.231. In the combined cost, the stem term (hundreds of units) swamps it. **The a-priori weights are therefore not neutral.** Proposed fix before Phase 4: normalise each soft term by its mean on whiten+RS oligos, so each weight is in "units of the unconstrained baseline". Then report a small weight grid or Pareto set rather than a single tuned point. **Needs a decision.**
  6. **Proxy check on real outputs** (2608 oligos): tier 2 (W = 80) vs full accessibility gives ρ = 0.95 pooled (0.79–0.999 per arm). Tier 1 gives ρ = −0.25, so it stays weak. Global MFE gives ρ = 0.32. The beam does not appear to exploit tier-2 blind spots.
  7. Encode cost: about 17–22 s per 4 KB for steering arms on 8 workers (tier-2 rescoring dominates), vs 0.1 s for baselines.
- **Steering codec:** implemented for alphabets B and A. Round trips are byte-identical. On 5 KB random data at 200 nt, B with P = 4–8 has no hard violations and needs ≤ 1 reroll. A at P = 8 needs a median of about 2 tries and up to 27, with one unrecoverable oligo, as F4 predicted. Encode time is about 0.05–0.12 s per oligo for B on one core (beam width 16).

## Part 0c — Revision 3 (Phase 2 addendum, binding)

### 0c.1 Arms
- **Co-primary steering arms: A** (quaternary + whitening + seeds) **and C** (run-length-limited quaternary, max run ≤ 3, plus steering). **B is demoted** to a comparison arm. Its role is to quantify the cost of run-length-1 enforcement (§0b.1).
- **Alphabet C** is a baseline construction, not a novelty claim. Its info positions use run-state mixed-radix coding: each info position carries a base-4 digit, except right after a run of `maxrun` equal bases, where it carries a base-3 digit over the bases ≠ previous. The decoder reads each radix from the preceding bases, including steering bases, so decoding is always possible, and max run ≤ 3 holds by construction.
  - Digits are extracted least-significant first from the whitened frame. The beam prunes any state whose remaining value cannot fit the remaining info positions, and requires remaining value = 0 at the end.
  - Frame size = 2·n_info − margin. The margin covers the expected n/16 forced positions plus 4 SD, at (2 − log₂3) bits each. That gives about 1.94 bits/nt at P = ∞, against the RLL(3) quaternary capacity of **1.98235 bits/nt** (largest eigenvalue 3.95137 of the run-state graph, computed).
  - An overflow (remaining value not 0) is a codec violation and triggers a seed reroll.
  - Prior RLL/GC-constrained DNA codes: Song, Cai, Zhang & Yuen 2018, *IEEE Commun. Lett.* 22(10):2004–2007, doi:10.1109/LCOMM.2018.2866566; Immink & Cai, arXiv:1812.06798; Dubé, Song & Cai 2019; Nguyen et al., *IEEE TIT* 2021 (see RELATED_WORK.md).
- **Seed exhaustion** (all seeds fail) is reported per arm. If A exhausts seeds at high P or on larger files, test a **4-nt seed field (81 seeds)** variant via `seed_trits: 4` (stored in the header; header oligos always use 3).
- **Phase 4 focus:** the high-density regime: A at P = 8–16, and C at matching P.

### 0c.2 Soft-cost normalisation, weight grid, ablations
- Each soft term (window-GC deviation, upper-tail R10 confusability, tier-1 structure proxy) is divided by its **IQR** of per-oligo totals on unconstrained whitened quaternary oligos (n = 400, fixed seed, scored by the beam's own scorer with hard constraints off; std recorded too). A weight of 1 therefore means "one baseline IQR". Implemented as `normalize: true` (`calibrate_spreads` in `codecs/steering.py`).
- **Weight grid (Phase 4, kept small because it multiplies Dorado runs):** w_conf × w_struct ∈ {0.5, 1, 2}², with w_gc = 1.
- **Single-term ablations:** gc = 0; conf = 0; struct = 0 (tier 1 off and tier-2 rescoring off).
- **Confusability is optimised against R10.4.1 only.** R9.4.1 is reported as a secondary metric.

### 0c.3 Bit-level packing
RS still runs over GF(256) bytes, on *rows* of ⌈d/8⌉ bytes, where d is an oligo's data bits. All rows of all blocks are concatenated row-major into one bitstream, and oligo q carries bits [q·d, (q+1)·d). Codecs now expose `frame_bits` instead of `frame_bytes`, and the CRC-16 is computed over the d data bits.
- **Effect:** density is strictly monotone in P, with no whole-byte plateaus (tested).
- **Cost:** a lost oligo erases one row in most columns, plus at most two extra bytes where it straddles a row boundary. Decoding groups columns by erasure pattern. Single undetected-bad-oligo location still works (tested).
- Header format version 2: `payload_bits` (u16) and `seed_trits` were added.

### 0c.4 Synthesis errors (addendum item d)
Sources checked are in `docs/notes/synthesis_errors.md`: Heckel et al. 2019, doi:10.1038/s41598-019-45832-6; Gimpel et al. 2023 and 2024; Lietard et al. 2021; Yeom et al. 2023; SOLQC. **No published, quantitative, synthesis-specific dependence of error rates on homopolymer length or GC content in our regime (runs ≤ 3, GC 40–60%) was found.** Heckel's homopolymer statement concerns *sequencing*. Only context-free per-nt rates exist; `dnastore/synthesis.py` implements them as sourced presets. These add the same error load to every arm, so they cannot reveal a synthesis-side benefit of constraints. **All conclusions about constraints vs. ECC are scoped to sequencing + amplification.**

### 0c.5 Coverage model reporting rule (addendum item e)
Every result that uses the PCR coverage model is reported twice: with the ASSUMED accessibility knob on, and with it set to 0 (`ChannelParams.zero_accessibility = True`).

### 0c.6 Accessibility reporting (decision 3)
Accessibility is a soft objective. Reported: the fraction of oligos with min anchor accessibility ≥ a* for a* ∈ {0.6, 0.686, 0.8}, and **between-arm comparisons by ranks and differences**. Every table states that absolute values depend on the ion model (ViennaRNA 2.7.2: monovalent Na⁺ correction only, no Mg²⁺).

### 0c.7 Phase 2 characterization sweep (single-payload, descriptive)
The sweep uses 4 KB of random data (seed 0), 200-nt oligos and the current normalized
steering scorer. These are **characterization results, not the equal-total-nt Phase 4
benchmark**; each arm has one payload and no channel replicates. Accessibility is at 60 °C
under the stated ViennaRNA monovalent-salt model. `conf tail` is excess R10 δ = 0.05 upper-
tail score; encode time is seconds on the configured 8 workers.

| Arm | Code density (bits/nt) | Violating oligos | Mean seed tries | Median min-anchor access | ≥0.686 | Conf tail | Encode (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| naive2bit | 1.882 | 87.8% | 1.01 | 0.766 | 69.1% | 4.79 | 0.58 |
| naive2bit-whiten | 1.882 | 86.3% | 1.02 | 0.709 | 59.0% | 4.74 | 0.54 |
| goldman | 1.463 | 3.4% | 1.00 | 0.764 | 64.6% | 4.42 | 0.57 |
| steering-A P=4 | 1.484 | 0% | 1.55 | 0.840 | 85.2% | 2.97 | 95.91 |
| steering-A P=6 | 1.605 | 0% | 2.39 | 0.839 | 83.3% | 3.31 | 48.57 |
| steering-A P=8 | 1.668 | 0% | 2.67 | 0.801 | 73.9% | 3.72 | 52.74 |
| steering-A P=12 | 1.733 | 0% | 3.65 | 0.814 | 76.0% | 4.16 | 58.50 |
| steering-A P=16 | 1.775 | 0% | 3.75 | 0.817 | 76.9% | 4.28 | 56.04 |
| steering-A P=∞ | 1.882 | 2.9% | 8.01 | 0.777 | 67.6% | 4.90 | 52.86 |
| steering-C P=4 | 1.435 | 0% | 1.01 | 0.864 | 81.9% | 2.72 | 51.65 |
| steering-C P=6 | 1.548 | 0% | 1.01 | 0.857 | 83.9% | 3.20 | 43.58 |
| steering-C P=8 | 1.619 | 0% | 1.02 | 0.842 | 78.3% | 3.62 | 42.96 |
| steering-C P=12 | 1.675 | 0% | 1.06 | 0.839 | 79.5% | 3.99 | 38.43 |
| steering-C P=16 | 1.699 | 0% | 1.06 | 0.839 | 78.6% | 4.12 | 38.29 |
| steering-C P=∞ | 1.812 | 0% | 1.37 | 0.805 | 74.3% | 4.80 | 30.47 |
| steering-B P=4 | 1.152 | 0% | 1.00 | 0.863 | 82.3% | 2.16 | 46.00 |
| steering-B P=6 | 1.246 | 0% | 1.00 | 0.854 | 83.3% | 2.61 | 43.89 |
| steering-B P=8 | 1.300 | 0% | 1.00 | 0.845 | 79.6% | 2.90 | 42.96 |
| steering-B P=12 | 1.359 | 0% | 1.00 | 0.842 | 83.2% | 3.16 | 41.52 |
| steering-B P=16 | 1.380 | 0% | 1.00 | 0.830 | 76.1% | 3.42 | 40.90 |
| steering-B P=∞ | 1.463 | 0% | 1.02 | 0.748 | 70.8% | 4.39 | 31.08 |

The table supports A and C as co-primary arms: A is denser at the same P, while C avoids
seed rerolls and has better anchor accessibility in this one-payload sweep. B is much less
dense, but its steering recovers accessibility relative to the run-length-1 Goldman baseline.
Steering lowers the reported R10 confusability tail at smaller P; the full Phase 4 weight
grid and channel results are still needed before judging whether that pays for its density
cost. The table and per-oligo source data are in `results/phase2/sweep_v2/summary.csv` and
`per_oligo.csv`.

### 0c.8 Phase 3 channel calibration and evaluation plan
Calibration against public R10.4.1 storage reads found 1.51% per-base error in real Dorado
hac, 1.55% in identity-matched Badread, and 9.12% in squigulator → Dorado hac; squigulator
also has a much larger short-run systematic floor. The calibration dataset is plasmid-fragment
data, not 200-nt oligos. Accordingly, **Badread is the main read-level evaluation channel**,
uniform IDS is the controlled second channel, and squigulator → Dorado is retained as a
miscalibrated stress test. Badread does not reproduce the small observed read-to-read
systematic component at homopolymer lengths ≥ 5. Full measurements and limits are in
`docs/notes/channel_calibration.md`.

Phase 3 updates the format to version 3 for CRC-32 and adds CRC-guided bounded repair plus
outer RS errors-and-erasures. The codec layout and common decoder are shared across arms;
results must show both `erasure` and `repair` strengths. This supersedes the version-2
format note in §0c.3.

### 0c.9 Phase 3 clustering scale check (completed 2026-10-08)

The 100,000-byte, six-arm run completed all 12 arm/channel cases at requested mean
depth 10. Pool sizes range from 3,678 to 5,742 oligos. All pairs of distinct bodies
were checked: minimum normalized edit distances are 0.38125 (whitened quaternary),
0.325 (Goldman), 0.375 (steering A), 0.36875 (steering C), 0.31875 (steering B) and
0.3625 (fountain); no pair is within 0.30L. Measured cluster merge rates and
`groups_multi_index` are zero throughout. Badread has zero split rate and all arms
recover at both decoder strengths. IDS9 has split rates of 2.04–6.41%; no arm recovers
the file at either strength. IDS9 here uses 4% substitutions, 2.5% insertions and 2.5%
deletions. These are one-payload clustering checks with unmatched synthesis budgets,
not the Phase 4 comparison. Source: `results/phase3/cluster_scale_100KB/cluster_scale.csv`.

## Part 1 — Flaws, ambiguities, and risks (blunt)

### F1. The confusability-weighted channel is circular. Results from it cannot support a readout claim.
If the encoder minimises metric *C* and the channel injects errors in proportion to *C*, the steering codec wins by construction. A reviewer will call this out immediately.
**Resolution:**
- (a) Report confusability-channel results only as a *sanity check* ("the optimiser does what it is told"). Never present them as evidence of better nanopore readability.
- (b) Make the real evidence an **independent** path: squigulator raw-signal simulation followed by an actual basecaller. Environment check: squigulator is **not on PyPI**, but prebuilt binaries are on GitHub and there is a bioconda recipe. There is an **RTX 3050 GPU**, so Dorado basecalling is plausibly feasible (4 GB VRAM may require the `fast`/`hac` models, not `sup`). I'll confirm in Phase 3 and say plainly if it fails.
- (c) Even that path is partly circular: squigulator generates signal from ONT's own k-mer level model, which is the table we optimise against. Disclose this. Wet-lab nanopore data is the only non-circular test.

### F2. The R10.4.1 table has no pA values and no standard deviations. "δ pA derived from the level std" is impossible for it.
I checked the repository (`nanoporetech/kmer_models`, commit `4e56dae`, 2023-06-06):
- `dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt`: 262,144 rows of `kmer<TAB>level`. Levels are in **standardised units** (mean ≈ 0, sd ≈ 1.0, range −2.0 to 2.96). There is no per-k-mer std and no pA scale.
- `legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model`: 6-mers with `level_mean` (pA), a per-k-mer `level_stdv` (≈1.2–1.5+ pA), and noise parameters.

**[DECISION D1] Which model, and what δ?** Default: **R10.4.1 9-mer as primary**, because it is current chemistry, with **δ in standardised units**, configurable, and default δ = 0.05. Measured on 3k random 9-mers × 27 Hamming-1 neighbours: the 5/10/25/50th percentiles of |Δlevel| are 0.009/0.018/0.051/0.127, so δ = 0.05 marks about the closest quarter of neighbours as confusable. This was fixed before seeing any codec results; a sensitivity sweep over δ will be reported. Run **R9.4.1 6-mer as a secondary model** where δ = *c*·level_stdv in pA (default *c* = 1) can be defined exactly as you specified. Record file, path, commit SHA, and SHA-256 in `results/provenance.json`.

### F3. Novelty is thin. See RELATED_WORK.md.
Seed retry is DNA Fountain screening / guided scrambling. Alphabet B is Goldman's code. Current-level-aware constrained coding was published in 2025 (Whritenour et al.). Weindel et al. argue the whole constrained-coding premise is inefficient. Consequences for the design:
- Add a **matched-budget baseline**: whitening + seeds, P = ∞, with the nucleotides saved by not steering spent on extra RS parity. If steering-A loses to this at equal total nt under the uniform channel, say so in RESULTS.md.
- Add Whritenour's **adjacent-k-mer contrast** as a second readout metric, measured on every codec.

### F4. Steering does almost nothing for homopolymers in alphabet A. The seeds do the work.
Quick Monte Carlo (20k trials, 150-nt whitened payload, homopolymer ≤ 3 only; ideal steering that never extends a run across its position):

| P | P(pass, one seed) | P(all 16 seeds fail) |
|---|---|---|
| 3 | 1.000 | 0 |
| 4 | 0.555 | 2.4e-6 |
| 8 | 0.301 | 3.3e-3 |
| 16 | 0.224 | 1.7e-2 |
| ∞ | 0.171 | 5.0e-2 |

A steering base can only break runs that *cross* its position. Runs inside a P-length info segment are fixed by the data. So **only P ≤ 3 guarantees the homopolymer constraint in A**, and P = 16 is barely better than no steering. Adding GC-window, MFE, and primer constraints will lower the pass rates further. Implications:
- The seed field needs ≥ 4 bits (16 seeds), and failure rates will be visible at P ≥ 8. They will be reported, not hidden.
- Steering's real job in A is GC, MFE, and confusability, not homopolymers. The paper should say this.
- Alphabet B guarantees run length ≤ 2 automatically (info base ≠ previous; a steering base can repeat at most one neighbour). B's steering is entirely about GC, MFE, and confusability.

### F5. In alphabet B, a steering base remaps the entire next segment. This is a feature, but it changes the search.
Rotating ternary: base_i = f(base_{i−1}, trit_i). Changing the steering base changes the first info base after it, which changes the next one, and so on through the segment. So each steering choice selects one of 4 *entirely different* realisations of the next P bases. That gives much more steering power than in A, where one choice changes one base. It acts like a local scrambler. Two consequences: the beam has to re-score whole segments per branch (cheap), and it is a genuinely interesting difference between A and B worth analysing in the paper.
Error propagation: one substitution in B corrupts two trits (the base itself and the next one's reference). This is the same as Goldman. The per-oligo CRC catches it.

### F6. Header bootstrapping: the header cannot be encoded with the codec it describes.
The header carries codec ID, P, and alphabet, but the decoder needs those to parse the header oligo. Separately, **×5 identical copies of index 0 are pointless**: identical sequences are one molecular species. They are synthesised and amplified as one population, have the same PCR bias, and drop out together. Five copies of the same string add no independent robustness beyond synthesis copy number.
**Resolution:** the header uses a **fixed, codec-independent format**: rotating ternary with a fixed P_hdr = 8 steering, fixed primers, and its own CRC32. It is written as **R distinct oligos** (indices 0…R−1 in a reserved header index space, each with different whitening or rotation offset) so their sequences differ and fail independently. The decoder takes any one valid copy, or majority-votes fields across copies.

### F7. A seed field "itself constraint-encoded" cannot be whitened. The same applies to the index.
The decoder must read the seed and index *before* it knows the keystream. Both are therefore encoded in a fixed **rotating-ternary sub-field with steering** (same as the header format) and covered by the per-oligo CRC. The keystream is `PRNG(global_seed, oligo_index, seed_idx)`, so identical data blocks do not produce identical oligos.

### F8. The orientation marker is redundant.
Asymmetric primers already identify orientation: a read starting with F is forward, and a read starting with R (i.e. ending in RC(F)) is reverse. A marker costs nt and adds nothing when primers are present. **Proposal: drop it.** Use primer matching (edlib, infix mode, with an error threshold) for orientation. Reads matching neither primer are counted and discarded. **[DECISION D2]** Default: drop.
Clarifying the layout convention: the physical oligo is `F + body + RC(R)`.

### F9. The MFE threshold and conditions are unspecified, and MFE is slow.
- Vienna 2.7.2 is installed. `RNA.params_load_DNA_Mathews2004()` works (verified). Defaults are 37 °C and 1 M Na⁺, which are not PCR or sequencing conditions. **Default:** 37 °C, recorded in config. Threshold set from the *empirical distribution of whitened random oligos with our primers* (e.g. its 10th percentile), plus a sweep. This avoids using a magic number from a paper with different lengths and primers.
- DNA MFE is **not** strand-symmetric (G·T wobble pairs become C·A on the reverse complement, and dangles differ). Both strands are folded, and the score is min(MFE_fwd, MFE_rc).
- **Measured: ~85 ms per fold for 200 nt** (single core). 1 MB → ~6k oligos × 2 strands × K rescoring candidates. With K = 8, that is ~96k folds ≈ 2.3 core-hours ≈ 17 min on 8 cores, *per configuration*. The P-sweep × 2 alphabets × 6 files is roughly 100+ configurations. **So:** in-beam I use a cheap structure proxy (count of reverse-complement k-mer pairs, k = 4–6, within the oligo plus primers). Only the final K beam survivors get real MFE. Benchmark files are capped at **≤ 100 KB** for the full sweep, with one ~1 MB run per codec at the chosen P for timing. **[DECISION D3]** OK to cap the sweep at 100 KB?

### F10. Hard-constraint failure policy is undefined.
"Never silently relax" still leaves a choice of what goes into the pool for a failed slot. Options: (a) emit the best-effort oligo, flag it, and count it; (b) emit nothing and let the outer RS treat it as a pre-erased column; (c) abort the encode.
**[DECISION D4]** Default: **(a) in benchmark mode** (so we can measure whether violating oligos actually decode worse), **(c) in production `encode` mode** unless `--allow-violations` is passed. Every violation is logged as (oligo index, which constraint, value).

### F11. Baseline fairness: ECC and layout must be shared, or the comparison is meaningless.
All codecs except fountain use **identical** primers, header scheme, index encoding, inner CRC, and outer RS parameters. Only the payload mapping differs. Fountain replaces the outer RS with LT redundancy and keeps the inner CRC. Its overhead is set so that the total nt matches the other codecs at the same target, not tuned separately.
**[DECISION D5] Goldman fidelity.** The original uses 4× overlapping segments with alternate RC, which is *redundancy*, not just a code. Default: implement `goldman` as **rotating ternary + our shared ECC stack** (fair code-vs-code comparison), and optionally `goldman_faithful` (4× overlap) as a historical reference point, clearly labelled.

### F12. Compression confounds density comparisons.
zstd applies to text and WAV, not JPEG, PNG, MP3, or MP4. If density is reported against original bytes, "per file type" results mostly measure compressibility. **Two metrics:** *code density* = stored (post-compression) bits / nt, which is what compares codecs, and *file density* = original bits / nt, which is what users care about. Compress iff the format is not in the skip list **and** the ratio is < 0.98. Always log the trial ratio.

### F13. Outer RS block structure needs specifying.
reedsolo works over GF(256), so n ≤ 255 symbols per codeword. 1 MB / ~45 bytes per oligo ≈ 23k oligos, so the data is split into blocks of ≤ 255 oligos. Each RS codeword runs *across* oligos (one byte from each oligo in the block). Parity = ⌈r·k⌉ oligos per block, with r configurable (default 15%). Loss of > n−k oligos in one block loses that block. This is reported per block. (galois GF(2^16) would allow one large block, but is much slower. Left as an option, not the default.)
Inner check: **CRC-16 per oligo** (8 nt in A; ⌈16/log₂3⌉ = 11 trits in B). Per-oligo bytes are fixed per codec so that the RS matrix is rectangular.

### F14. Ternary ↔ binary conversion has to be specified, because log₂3 is asymptotic.
Per-oligo big-integer conversion: the payload's bits are treated as one integer and written in base 3 with a fixed trit count. Rate = ⌊L_info·log₂3⌋ bits per oligo, which is close to log₂3 with no cross-oligo propagation. Fixed-block alternatives (11 trits ↔ 17 bits = 1.545 bits/trit) lose 2.5%. Default: big-int.

### F15. Decoder with indels: "skip steering positions" only works after alignment.
Steering positions are fixed offsets, so a single indel shifts everything after it. Order of operations: primer trim → read the index field (fixed offset; errors here cause misgrouping, and the edlib fallback clusters by whole-body similarity) → **consensus first** (length-filtered majority vote over reads of the expected length; if too few, align reads to a medoid with edlib and vote per column) → then strip steering → decode → CRC. A consensus that fails CRC becomes an erasure. RS can also correct some undetected errors (2 errors cost as much as 1 erasure).

### F16. Statistics: 20 seeds per point is weak.
Binary success with n = 20 gives a 95% Wilson CI of about ±0.2 at p = 0.5. Default: **≥ 20 as the floor, 100 where affordable** (decode is cheap; encode is done once per configuration and reused across channel seeds). Report Wilson intervals. Also report **continuous** outcomes (fraction of oligos recovered before RS, residual erasures per block), which carry more information per run than pass/fail.

### F17. "dsDNA realities" vs ssDNA MFE.
Hairpins form when strands are single-stranded: during synthesis, at PCR denaturation/annealing, and in nanopore translocation. MFE at 37 °C is a proxy for all of these. That is fine, but say it is a proxy. Primer cross-talk is checked for F, R, RC(F), and RC(R) against all body windows with ≤ k mismatches (default k = 3 for 20-nt primers), plus 3′-end-anchored partial matches (the last 8 nt, ≤ 1 mismatch), because 3′ complementarity drives mispriming. In steering codecs, primer avoidance is a **hard constraint inside the beam**. In baselines it is screened and reported.

### F18. Sample files.
Synthetic images (ffmpeg `testsrc`) compress unrealistically well and have odd byte statistics. Default: one **real public-domain file per type**, downloaded with URL, licence, and SHA-256 recorded (e.g. Wikimedia Commons PD images, a Project Gutenberg text, LibriVox PD audio). ffmpeg (present) is used to transcode, trim, and create MP3/MP4/WAV from PD sources. All files 10 KB–100 KB for the sweep.

---

## Part 2 — Proposed concrete design

### 2.1 Oligo layout (default 200 nt)
```
[F 20][ index+seed sub-field (ternary, steered) ][ payload body (codec-specific, incl. steering) ][ CRC-16 ][ RC(R) 20 ]
```
For steering-A/B, steering positions run through the index sub-field, body, and CRC with a single period P. Their positions depend only on (P, layout), and the layout is stored in the header. Pattern within a steered region: `[P info][1 steer][P info][1 steer]…`, so steered rate = P/(P+1) × (2 or log₂3) bits per position.

### 2.2 Steering beam search
- **State:** prefix so far, plus incremental counters (current run, GC count in the window, a running cost vector).
- **Branching:** 4 choices at each steering position. Info segments are deterministic given the state (A: the whitened bits; B: trits plus the previous base).
- **Hard constraints (prune):** run ≤ 3; window GC within limits; global GC is checked at the end (its feasibility is checked on the way via a bound on the remaining positions); no primer hit.
- **Soft cost:** w_gc·|GC − 0.5| + w_struct·(RC-kmer-pair proxy) + w_conf·(confusability, both strands).
- **Width:** B (default 16). The final K = min(B, 8) survivors are rescored with real MFE (fwd and RC, with primers). The best one that passes the MFE threshold is chosen. If none pass in A, retry with the next seed. If none pass in B (or A after all seeds), apply the policy in D4.

### 2.3 Confusability (exact definition, to be implemented in `screening.py`)
For model M with levels ℓ(x) for k-mers x ∈ Σ^k: conf(x) = |{ y : d_H(x,y) = 1, |ℓ(x) − ℓ(y)| < δ }| / (3k), which lies in [0, 1].
Oligo score = mean and max of conf over all k-mers of the oligo **and of its reverse complement** (both strands get read). The beam uses the mean. Both are reported.
Second metric (Whritenour-style): adjacent contrast = min and mean over i of |ℓ(x_i) − ℓ(x_{i+1})|.

### 2.4 Header (fixed format, see F6)
magic `DNAS` (4 B) · version (1 B) · file-type code (1 B; magic-byte detection) · original length (8 B) · compressed flag (1 B) · stored length (8 B) · codec ID (1 B) · codec params (P, alphabet, ECC k/n, block count, layout lengths; TLV, ≤ 24 B) · global PRNG seed (4 B) · CRC32 of original (4 B) · header CRC32 (4 B). Split over the R header oligos if needed (each contains all fields, so any one copy suffices).

### 2.5 Metrics definitions
- **code density (net)** = stored data bits / body nt (body includes steering; excludes primers, index/seed, CRC, header oligos, parity oligos).
- **effective density** = original file bits / all synthesised nt (everything).
- Also report stored-bits / all nt.
- **mass** = Σ oligos × copies × length × 330 g/mol / N_A (ssDNA), and 650 g/mol per bp for dsDNA. Copies come from config.

### 2.6 Reproducibility
One command: `dnastore benchmark configs/paper.yaml` → `results/<run_id>/` containing CSVs, figures (PNG + SVG), a `provenance.json` (git SHA, package versions, k-mer table SHA, sample file SHAs), and the config copy. Every RNG is derived from a single root seed via `numpy.random.SeedSequence`.

### 2.7 Environment notes (verified)
Python 3.13.11 (spec says 3.11+, OK). Installed: ViennaRNA 2.7.2, numpy 2.2.6, zstandard 0.24, typer 0.20. ffmpeg is present. GPU: RTX 3050 Laptop. Not yet installed: reedsolo, edlib, hypothesis, pandas/matplotlib (to check), squigulator. I'll create a dedicated venv/conda env with pinned versions in Phase 1.

---

## Decisions (resolved in Revision 1: defaults accepted)
| ID | Question | Default if you say "go" |
|---|---|---|
| D1 | k-mer model and δ units | R10.4.1 9-mer primary (δ = 0.05 std. units ≈ 25th pct of neighbour abs(Δlevel)); R9.4.1 6-mer secondary (δ = 1·level_stdv pA) |
| D2 | Drop the orientation marker? | Drop; use primer-based orientation |
| D3 | Cap sweep files at 100 KB? | Yes; one ~1 MB run per codec for timing |
| D4 | Failure policy | Benchmark: emit + flag + count; encode CLI: abort unless `--allow-violations` |
| D5 | Goldman fidelity | Code + shared ECC; optional `goldman_faithful` |
| — | Framing | Integration/trade-off study, with the matched-budget "embrace errors" baseline included |
