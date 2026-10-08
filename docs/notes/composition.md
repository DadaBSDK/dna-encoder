# Why rotating ternary folds more stably: composition analysis

Script: `scripts/composition_analysis.py` (seed 20261001). Outputs: `results/phase2/composition/`, containing `summary.json`, `mfe_summary.csv`, `mfe_per_sequence.csv`, `comp_kmer_pairs.csv`, `dinucleotide_freq.csv`, and `composition_mfe.{png,svg}`.
Setup: 160-nt sequences without primers, n = 2000 per source. ViennaRNA 2.7.2 with DNA Mathews 2004 parameters, default salt, at 37 °C and 60 °C. A smaller set of 300 per source was run with primers. Nearest-neighbour stacking energies come from SantaLucia (1998) unified dG37 values, *PNAS* 95:1460, doi:10.1073/pnas.95.4.1460. Complementary k-mer pairs are counted only between non-overlapping windows separated by a loop of at least 3 nt.

## Sources
| source | definition |
|---|---|
| `iid` | uniform i.i.d. quaternary (the whitened-payload model) |
| `iid_hp3` | `iid` rejection-filtered to max homopolymer ≤ 3 (the realistic constraint) |
| `p0.15`, `p0.10`, `p0.05` | Markov chain: P(next = prev) = r, other 3 bases uniform |
| `rot` | rotating ternary with uniform trits (r = 0) |
| `shiftban` | control: forbids a→next(a) (A→C, C→G, G→T, T→A) instead of a→a. Same transition collision probability as `rot` (1/3), homopolymers allowed. |

## Results (37 °C; 60 °C in brackets)
| source | median MFE (95% CI) | paired frac | stacks/seq | GC-stack frac | mean stack dG | comp. 6-mer pairs (obs / theory) |
|---|---|---|---|---|---|---|
| iid | −17.2 (−17.4, −17.0) [−8.8] | 0.52 | 29.9 | 0.33 | −1.51 | 2.64 / 2.62 |
| iid_hp3 | −17.3 (−17.5, −17.1) [−8.9] | 0.53 | 30.2 | 0.32 | −1.51 | 2.68 / — |
| p0.15 | −18.5 [−9.7] | 0.54 | 31.4 | 0.28 | −1.50 | 3.26 / 3.40 |
| p0.10 | −20.7 [−10.7] | 0.57 | 33.1 | 0.25 | −1.49 | 4.62 / 4.62 |
| p0.05 | −23.3 [−12.4] | 0.60 | 35.2 | 0.23 | −1.48 | 6.81 / 6.89 |
| **rot** | **−27.1 (−27.3, −26.9)** [−14.7] | 0.63 | 38.0 | 0.22 | −1.47 | 10.97 / 11.04 |
| shiftban | −25.0 (−25.2, −24.8) [−13.5] | 0.59 | 35.9 | 0.30 | −1.46 | 11.08 / 11.04 |

With primers (200 nt, 37 °C): iid −22.85 (−23.4, −22.4), rot −31.8 (−32.35, −31.2).

Dinucleotides: theory and observation agree to ±0.001. For `rot`, XX = 0 and XY = 1/12. Self-complementary steps (AT, TA, CG, GC) make up 0.333 of `rot` vs 0.250 of `iid`. Strong steps (CG, GC, GG, CC) make up 0.167 vs 0.250. The mean NN dG per step in the sequence is almost the same (−1.401 vs −1.406 kcal/mol).

Within each source, the per-sequence Pearson r between complementary 6-mer pairs and MFE is −0.35 to −0.43.

## Mechanism
- **(i) More complementary k-mers: supported, and it is the main driver.** For any reverse-complement-symmetric first-order chain, the chance that two distant k-windows are reverse complements is ¼·c^(k−1), where c = Σ_b T(a,b)² is the transition collision probability. Forbidding repeats raises c from 1/4 to 1/3. That gives (4/3)^(k−1) more hairpin-capable complementary pairs: 4.2× at k = 6 and 7.5× at k = 8. The observed counts match theory exactly. MFE falls monotonically as the repeat probability drops (−17.2 → −18.5 → −20.7 → −23.3 → −27.1).
- **(ii) More strong (GC) stacks: rejected.** `rot` has *fewer* GG/CC/CG/GC steps, because GG and CC are banned. Its MFE structures have a lower GC-stack fraction (0.22 vs 0.33) and weaker average stacks (−1.47 vs −1.51 kcal/mol). The extra stability comes from *more* stacked pairs (38 vs 30 per sequence; 63% vs 52% of bases paired), not stronger ones.
- **Control.** `shiftban` has the same collision probability but *allows* homopolymers. It reproduces about 79% of the gap: (25.0 − 17.2)/(27.1 − 17.2). So the cause is mostly the **entropy reduction of the first-order chain**, not the absence of repeats as such. The remaining ~2.1 kcal/mol fits `rot`'s enrichment in self-complementary steps (0.333 vs 0.167 in `shiftban`). That attribution is plausible but not separately tested.
- **≤ 3 vs ≤ 1.** Filtering i.i.d. sequences to max homopolymer ≤ 3 costs **nothing measurable**: −17.25 vs −17.2, with overlapping CIs, and the same complementary-pair count. The cost belongs to codes that remove repeats at *every* position (run ≤ 1, as in rotating ternary), not to the usual RLL-3 constraint.

## Candidate result (for DESIGN.md)
> **Homopolymer constraints can carry a hidden structure cost, and the size depends on how the constraint is enforced.** Rotating ternary (run length 1) folds about 10 kcal/mol more stably than whitened quaternary: median MFE −27.1 vs −17.2 for 160 nt at 37 °C, −14.7 vs −8.8 at 60 °C, and −31.8 vs −22.9 for full 200-nt oligos with primers. The cause is lower entropy in the first-order chain. Banning one of four successors raises the transition collision probability from 1/4 to 1/3, which multiplies the expected number of hairpin-capable reverse-complementary k-mer pairs by (4/3)^(k−1): 4.2× for 6-mers, and theory matches observation. GC stacking is not the cause; `rot` has fewer GC stacks. A control that bans a different successor reproduces ~79% of the gap. Rejection-filtering to the conventional max run ≤ 3 has no measurable cost. **Implication:** alphabet B starts from a structural deficit that its steering bases must pay back. Any comparison of B against whitened quaternary has to charge B for it.

## Caveats
- These are MFE values only: no ensemble or accessibility measures, which are the primary metric from Phase 2 onward. Default ViennaRNA salt.
- Uniform trits stand in for whitened data. Sequences without primers are 160 nt; the with-primer set is smaller (n = 300).
- The ~21% residual attributed to self-complementary step enrichment is an inference.
- Dangles and the ViennaRNA DNA parameter set may change absolute values. Only the comparisons are claimed.
