# Channel simulator and noisy-read primitives (Phase 3)

Modules: `dnastore/channel.py`, `dnastore/consensus.py`. Tests: `tests/test_channel.py`, `tests/test_consensus.py`.

## Channel model (`simulate(oligos, ChannelParams, rng, features)`)
1. **Copy weights `w_i`:**
   - `poisson`: all equal.
   - `lognormal`: LogNormal(−σ²/2, σ), mean 1.
   - `model`: `dnastore.coverage.expected_copies(features, SENSITIVITY_GRID[preset])`. `features["gc"]` is computed from the oligos if absent; `features["accessibility"]` is optional.
   - **`zero_accessibility=True` sets the ASSUMED accessibility penalty to 0 and drops the accessibility feature.** Every coverage-model result must be reported with this switch both off and on.
2. **Dropout:** each oligo is lost independently with probability `dropout`.
3. **Depth:** the total read count is fixed at N = round(`mean_coverage` · n_oligos) and allocated by Multinomial(N, w/Σw). This models a sequencing run of fixed depth, so dropout redistributes reads to survivors and does not reduce N.
4. **Orientation:** each read is reverse-complemented with probability `rc_frac` (default 0.5). **Errors are applied in read orientation**, to the strand as sequenced.
5. **IDS errors,** independently at each template position j: with u ~ U(0,1), delete if u < p_del; else substitute (uniform over the 3 other bases) if u < p_del + p_sub. Independently, insert one uniform random base before j with probability p_ins. Per-read counts (`n_sub`, `n_ins`, `n_del`) are recorded, and len(read) = L + n_ins − n_del.

## Confusability-weighted mode (`conf_weighted=True`)
> **PROXY for nanopore behaviour, NOT a validated error model, and CIRCULAR with respect to the steering encoder's confusability objective** (the encoder minimises the same table this mode uses to place errors). Results under this mode are a **sanity check only** and must be labelled circular in every table and figure.

- Rates become p_x,j = p_x · m_j, with m_j = c̄_j / μ.
  - c̄_j is the mean confusability (R10.4.1 9-mer, δ = 0.05 by default) of the k-mers lying fully inside the read-orientation strand that cover position j.
  - μ is the mean of the per-k-mer table.
- **Normalisation is exact in expectation.** For a uniform random sequence each covering k-mer is uniform over Σ^k, so E[c̄_j] = μ and E[m_j] = 1 at every position, ends included. Test result: empirical substitution rate 0.02 ± 0.002 at base 0.02, mean multiplier 1.00 ± 0.02.
- Clipping: if p_del + p_sub > 0.9 both are scaled down proportionally; p_ins ≤ 0.9. Clipping never binds at realistic rates; m_j ranges about 0.5–2.
- Multipliers are computed on revcomp(oligo) for RC reads.

## Read processing (`consensus.py`)
- **`locate_and_orient(read, F, tail, max_primer_ed=4, margin=80)`:** edlib HW (infix) search for F in the first len(F) + 80 nt and for `tail` = rc(R) in the last len(tail) + 80 nt, on both the read and its reverse complement. The expanded window supports adapter-flanked reads. Both primers must match within 4 edits, and the orientation with the lower total edit distance wins. Returns `(body, is_rev, {ed_forward, ed_tail})`. The body length varies with indels.
- **`consensus(bodies, expected_len)`:**
  - Every read is aligned (edlib NW path) to a reference. The reference is the column majority of exact-length reads if there are ≥ 5 of them making up ≥ 30% of the reads, otherwise a medoid of ≤ 15 sampled reads.
  - Votes are taken per reference column (a base or "−" for deletion) and per gap (the inserted string, or none). The process iterates for ≤ 2 rounds.
  - The output is not padded or trimmed to `expected_len`; a wrong length normally fails CRC downstream.
  - *Column majority alone was rejected.* At 3% IDS, most exact-length reads still contain a compensating insertion/deletion pair, and the shifted segment corrupts the column vote. Exact recovery fell from 0.965 at 10× to 0.735 at 20×. A regression test checks this.
- **`cluster(bodies, max_ed_frac=0.15)`:** greedy assignment, in random order, to the first centre within ⌊0.15 · len⌋ edits (edlib with a k bound); otherwise the read starts a new centre. This is the fallback for reads whose index cannot be decoded.

## Validation (160-nt bodies, IDS split equally across sub/ins/del, 200 trials per cell, rc_frac = 0)
Exact-recovery fraction of the consensus (mean edit distance per nt in brackets):

| total IDS rate | 3× | 5× | 10× | 20× |
|---|---|---|---|---|
| 1% | 0.945 (3e-4) | 0.995 (0) | 1.000 (0) | 1.000 (0) |
| 3% | 0.570 (3.6e-3) | 0.930 (4e-4) | 0.995 (0) | 1.000 (0) |
| 6% | 0.175 (1.3e-2) | 0.545 (3.7e-3) | 0.935 (5e-4) | 1.000 (0) |
| 10% | 0.000 (3.4e-2) | 0.110 (1.4e-2) | 0.560 (3.9e-3) | 0.900 (7e-4) |

Coverage here is the *mean* (multinomial sampling), so some trials get fewer reads. Consensus costs about 10–65 ms per group.

## Throughput (one core)
- `simulate`: 1e5 reads of 200 nt in 7.9 s; 1e6 in 49.7 s; conf-weighted 1e5 in 27 s (multiplier precomputation is per oligo).
- `locate_and_orient`: 20k reads in 0.9 s. Primers found in 99.96% at 3% IDS, orientation correct in 100%.

## Assumptions and limitations
- Errors are i.i.d. per position (except in conf-weighted mode). There are no homopolymer-length errors, burst errors, or position-dependent rates beyond conf weighting, and no chimeras or synthesis errors. Synthesis and PCR substitutions are out of scope unless a published model is added.
- The total depth is fixed, not Poisson.
- Calibrated Badread is the main read-level evaluation; uniform IDS is the controlled reference. Squigulator + Dorado is a miscalibrated raw-signal stress test (see `channel_calibration.md` and DESIGN.md §0c.8).
