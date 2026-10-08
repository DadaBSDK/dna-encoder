# LT fountain core (`dnastore/lt.py`): notes

Our own implementation, written from Luby (2002), "LT codes" (FOCS 2002, doi:10.1109/SFCS.2002.1181950), and MacKay, *Information Theory, Inference, and Learning Algorithms*, ch. 50. **No code from the DNA Fountain repository was consulted.** The only thing reused from Erlich & Zielinski (2017, doi:10.1126/science.aaj2038) is their published robust-soliton parameter *choice* (c = 0.025, δ = 0.001), which is our default and is cited as theirs.

## Algorithm
- **Degree distribution:** robust soliton (MacKay eqs. 50.4–50.6).
  - ρ(1) = 1/k, ρ(d) = 1/(d(d−1)).
  - S = c·ln(k/δ)·√k.
  - τ(d) = S/(kd) for d < k/S, and τ(k/S) = S·ln(S/δ)/k (spike index rounded to the nearest integer in [1, k]).
  - μ = (ρ + τ)/Z.
  - Sampling uses an inverse-CDF lookup (`searchsorted`).
  - Mean degree at the default c and δ: k = 100 → 10.5, 1000 → 15.4, 10 000 → 20.0, 25 000 → 21.8.
- **Droplet = (32-bit seed, payload).** The seed seeds numpy PCG64, which draws the degree, then d distinct neighbours with Floyd's algorithm (O(d), independent of k). The payload is the XOR of the neighbour segments. The decoder regenerates neighbours from the seed, so it is deterministic for a pinned numpy version.
- **Decoding:**
  1. Peeling (belief propagation with a ripple of degree-1 droplets).
  2. If peeling stalls, Gauss–Jordan elimination over GF(2) on the residual system: bit-packed coefficient rows, with payload rows XORed in vectorised numpy. Only the unknowns with a fully reduced pivot row are declared solved.
  - Duplicate seeds are ignored and counted.
  - Failure is reported explicitly (`ok=False`, `n_missing`, `segments=None`), never returned as partial data.
- **Screening hook:** `generate_screened(segments, seed_iter, accept, n_needed, dist, stats)` draws seeds in order and keeps droplets passing `accept(seed, payload)`. It records drawn, accepted and rejected counts, and raises if the seed iterator is exhausted.

## API
```python
RobustSoliton(k, c=0.025, delta=0.001)         # .pmf .cdf .sample(rng) .mean_degree .S .Z
droplet_neighbors(seed, k, dist) -> list[int]  # sorted, deterministic
encode_droplet(segments, seed, dist) -> bytes
generate_screened(segments, seed_iter, accept, n_needed, dist=None, stats=None) -> list[(seed, bytes)]
LTDecoder(k, seg_len, dist=None, max_ge_unknowns=20000)
    .add(seed, data); .decode() -> LTResult(ok, segments, n_missing, n_peeled, n_ge,
                                            n_droplets, n_duplicates, info)
```

## Measured overhead vs. success
Setup: 50 trials per cell, distinct random seeds per trial (`rng.choice(2**32, n, replace=False)`), 8-byte segments (success depends only on the neighbour graph, not on segment length), default c and δ. Overhead ε means n = round(k(1+ε)) droplets received with no losses.

| k | ε = 0 | 2.5% | 5% | 7.5% | 10% | 15% | 20% |
|---|---|---|---|---|---|---|---|
| 100 | 8/50 | 24/50 | 48/50 | 50/50 | 50/50 | 50/50 | 50/50 |
| 1 000 | 7/50 | 50/50 | 50/50 | 50/50 | 50/50 | 50/50 | 50/50 |
| 10 000 | 3/50 | 50/50 | 50/50 | 50/50 | 50/50 | 50/50 | 50/50 |

Number of trials (of 50) in which the GE fallback solved at least one segment:

| k | ε = 0 | 2.5% | 5% | 7.5% | 10% | 15% | 20% |
|---|---|---|---|---|---|---|---|
| 100 | 50 | 50 | 50 | 50 | 50 | 50 | 50 |
| 1 000 | 50 | 50 | 50 | 50 | 50 | 45 | 17 |
| 10 000 | 50 | 50 | 50 | 22 | 3 | 0 | 0 |

So at ε ≥ 2.5% (k ≥ 1000) or ε ≥ 7.5% (k = 100), every trial succeeded. **Much of that low-overhead success is due to the GE fallback, not to peeling alone.** Pure peeling with these parameters needs noticeably more overhead, especially at small k. This is a decoder property and should be stated whenever fountain results are compared with RS. Median missing segments at ε = 0 were 46 (k = 100), 550 (k = 1000) and 7099 (k = 10 000).

Timing (k = 25 000, 30-byte segments, 8% overhead, one core, other processes running): encoding 27 000 droplets ≈ 6 s, decoding ≈ 6 s (peeling only).

## Honest notes and pitfalls
- **Screening that is linear in the payload breaks decodability.** If `accept` is a GF(2)-linear function of the payload (for example, the parity of one bit), every accepted droplet's neighbour set has even overlap with the set of segments where that bit is 1. All accepted rows are then orthogonal to that indicator vector, so rank ≤ k−1 however many droplets are collected (`test_linear_screening_rule_makes_system_rank_deficient`). DNA screening rules (GC window, homopolymers, motifs on the mapped bases) are not linear, so they don't cause exact rank loss. They do bias which neighbour sets survive. The codec should report the screening rejection rate and **measure** decoding overhead under screening rather than assuming the unscreened curve above. Whitening the segments before LT encoding reduces data dependence of the accept decisions.
- The decoder assumes droplets are uncorrupted (the inner CRC is upstream). A corrupted droplet that passes CRC would propagate silently through peeling. The outer integrity check (file CRC32) catches it at the end, but it is not correctable here.
- GE cost grows as O(u² · m / 8) for u residual unknowns. `max_ge_unknowns` (default 20 000) caps it. If peeling leaves more unknowns than that, GE is skipped and the result is reported as a failure.
- Determinism relies on numpy's PCG64 and `Generator.random` / `integers` streams, which are stable for the pinned numpy 2.2.6.
- Seeds are 32-bit. For DNA, the codec stores the seed in the oligo, so the seed field length (16 nt quaternary in DNA Fountain; to be chosen in `codecs/fountain.py`) is a density cost to account for.
