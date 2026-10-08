# Related Work

Status: Phase 0 scan, done 2026-10-01 with web search. This is **not** a systematic review. I read abstracts, repository READMEs, and summaries. I read full text for only a few papers, and those are marked. Before any novelty claim goes into the paper, read the cited papers in full and run a proper search (Google Scholar, plus IEEE ISIT/TIT, *Nucleic Acids Research*, *Nature Communications*, and bioRxiv from 2019 to 2026).

---

## ⚠️ Overlap warning (read first)

**Little of this project's design is new taken one piece at a time.** The parts overlap with prior work as follows:

1. **The nanopore confusability objective overlaps directly with Whritenour, Civelek & Farnoud (2025).** They constrain DNA sequences so that the ONT current levels of **adjacent** k-mers differ by more than τ pA. They use the R9.4.1 6-mer table with a pruned de Bruijn graph and a state-splitting encoder. They report about 6× lower edit distance with a Viterbi decoder, and 24–25% fewer errors from Guppy/Dorado, all in simulation. This is the closest prior work, and it **must** be cited and compared against. What differs in our plan: we score **Hamming-1 neighbours** with similar levels (closer to substitution confusion than to segmentation), we use the R10.4.1 9-mer table, and we enforce the objective softly through sparse steering bases rather than a hard constrained code. Those differences are modest. A reviewer will ask why our metric beats theirs, so we should implement their adjacent-k-mer contrast as a second metric and report both.
2. **Steering bases are a known technique.** Inserting redundant symbols so that a sequence meets constraints is the core of Knuth-style balancing (Dubé, Song & Cai 2019), of guided scrambling, and of the greedy iterative coder by Park, Lee & No (2021). Alphabet A's mechanism (XOR whitening, then retrying seeds until screening passes, then storing the seed) is the screening loop of **DNA Fountain** (Erlich & Zielinski 2017) and is guided scrambling in all but name. Alphabet B is the **Goldman (2013)** rotating ternary code.
3. **MFE as an encoding objective has been done before.** Cao et al. (2021) built coding sets under an MFE constraint. Secondary-structure-avoidance codes (Nguyen et al.) handle structure combinatorially. The yin–yang codec (Ping et al. 2022) searches over segment pairings to pass GC and homopolymer screening (its reference code exposes `search_count`).
4. **A published result argues against the premise.** Weindel, Gimpel, Grass & Heckel (arXiv 2023 / IEEE TMBMC) argue that *randomising, accepting errors and adding ECC redundancy* is more efficient than constrained coding in realistic error regimes. If that holds, steering bases spend nucleotides that would do more good as parity. **Our benchmark must include this comparison directly**: whitening only (P = ∞) plus extra RS parity at the *same total nt* as steering-A/B.

**What might still be new (unverified; I did not find these combined):** (a) sparse, fixed-position steering bases chosen by beam search over a *joint* cost of GC, homopolymer, MFE of the full oligo including primers on *both strands*, and k-mer confusability; (b) how steering interacts with rotating ternary, where one steering base remaps the whole following segment, so it acts as a local scrambler rather than a single-symbol fix; (c) a Hamming-neighbour confusability score derived from the R10.4.1 9-mer level table. The honest framing is an **integration and empirical trade-off study** (density vs. structure vs. readout), not a new coding primitive.

---

## 1. Foundational systems

| Work | What it does | Overlap with this project |
|---|---|---|
| **Church, Gao & Kosuri (2012)**, *Science* 337:1628. [doi:10.1126/science.1226355](https://doi.org/10.1126/science.1226355) | 1 bit per base (A/C = 0, G/T = 1). The free choice is used to avoid homopolymers and balance GC. Addressed oligos; 5.27 Mbit book. | **Conceptually the closest ancestor of steering**: it gives up rate (1 bit/nt instead of 2) to gain freedom to satisfy constraints. Our steering is a sparse, optimised version of the same trade. |
| **Goldman et al. (2013)**, *Nature* 494:77–80. [doi:10.1038/nature11875](https://doi.org/10.1038/nature11875) | Huffman coding to trits, then a rotating ternary code (base ≠ previous base, so no homopolymers). 100-nt segments overlapping 4× with offset 25, alternate segments reverse-complemented, plus parity trit and index. | Alphabet B **is** this code. Our baseline `goldman.py` should state whether it reproduces the 4× overlap (see DESIGN.md §D3). |
| **Erlich & Zielinski (2017)**, "DNA Fountain", *Science* 355:950–954. [doi:10.1126/science.aaj2038](https://doi.org/10.1126/science.aaj2038) | Luby-transform droplets with an RS inner code. Each candidate droplet is screened (GC, homopolymer) and discarded if it fails, so the droplet seed acts as the free parameter. | Alphabet A's seed retry is the same screening idea. Our `fountain.py` is a reimplementation of the concept, not of their code. |

## 2. Constrained, guided-scrambling, and search-based codes

| Work | What it does | Overlap |
|---|---|---|
| **Immink & Cai**, "Properties and constructions of constrained codes for DNA-based data storage". [arXiv:1812.06798](https://arxiv.org/abs/1812.06798) | Capacity analysis and constructions for run-length and GC constraints. | Gives the capacity bounds for comparison: the maximum rate under RLL-3 + GC, against which our steering rate P/(P+1)·2 should be compared. |
| **Song, Cai, Zhang & Yuen (2018)**, "Codes with run-length and GC-content constraints for DNA-based data storage", *IEEE Commun. Lett.* 22(10):2004–2007. [doi:10.1109/LCOMM.2018.2866566](https://doi.org/10.1109/LCOMM.2018.2866566) | Efficient codes satisfying run-length and GC constraints together. | Our alphabet C (run-length-limited quaternary, max run ≤ 3) is a simple baseline from this family, not a contribution. RLL(3) quaternary capacity is 1.98235 bits/nt (computed). |
| **Dubé, Song & Cai (2019)**, "DNA codes with run-length limitation and Knuth-like balancing of the GC contents", SITA 2019. [Semantic Scholar](https://www.semanticscholar.org/paper/d070b0b431e3cfdd46fde6d356f6bb0932e42683) | Embeds data under RLL-3, then balances GC with a Knuth-style prefix and a small redundant index. | Redundant symbols spent to satisfy constraints, with a *guaranteed* outcome. Steering is heuristic (beam search) and therefore has no guarantee. Be ready to argue why. |
| **Nguyen, Cai, Immink & Kiah**, "Capacity-approaching constrained codes with error correction for DNA-based data storage", *IEEE Trans. Inf. Theory* (2021). [ResearchGate](https://www.researchgate.net/publication/350148599_Capacity-Approaching_Constrained_Codes_With_Error_Correction_for_DNA-Based_Data_Storage) | Constrained codes combined with ECC, near capacity. | Strong theory baseline for the GC + RLL part of our objective. |
| **Park, Lee & No (2021)**, "Iterative DNA coding scheme with GC balance and run-length constraints using a greedy algorithm". [arXiv:2103.03540](https://arxiv.org/abs/2103.03540) | Greedy iterative encoder; 1.8523 bits/nt for m = 3, α = 0.05. | Greedy search over free choices to meet GC + RLL. Our beam search is a wider version of the same idea with more objectives. **Compare densities directly**: 1.85 bits/nt is a high bar for steering-A at P ≥ 12. |
| **Ping et al. (2022)**, yin–yang codec, *Nat. Comput. Sci.* [doi:10.1038/s43588-022-00231-2](https://doi.org/10.1038/s43588-022-00231-2); [code](https://github.com/ntpz870817/DNA-storage-YYC) | Two rules map 2 bits to 1 nt. Pairs binary segments and searches (`search_count`) for pairings that pass GC (`max_content`) and homopolymer screening. I could not confirm the exact free-energy screening step from the README; the paywall blocked the full text. | Search over free parameters to meet constraints. Also tested with real files in 200-nt pools, which is a useful reference point. |
| **Löchel et al. (2022)**, "Fractal construction of constrained code words for DNA storage systems", *NAR* 50:e30. [doi:10.1093/nar/gkab1209](https://doi.org/10.1093/nar/gkab1209) | Chaos-game-representation construction of codebooks under GC, homopolymer, and motif constraints. | Codebook alternative to steering. |
| **Welzel et al. (2023)**, DNA-Aeon, *Nat. Commun.* 14:628. [doi:10.1038/s41467-023-36297-3](https://doi.org/10.1038/s41467-023-36297-3) | Arithmetic coding over a constrained model (GC, homopolymer, undesired motifs) concatenated with ECC; corrects substitutions, indels, and dropouts. | A modern, strong constrained-coding baseline. We do **not** benchmark against it in the current plan, and a reviewer may ask why. |
| **Press et al. (2020)**, HEDGES, *PNAS* 117:18489. [doi:10.1073/pnas.2004821117](https://doi.org/10.1073/pnas.2004821117) | Hash-based inner code that corrects indels directly and allows sequence constraints. | Our inner layer is CRC-only (detection). HEDGES shows that inner indel *correction* is a known and stronger alternative. |
| **Weindel, Gimpel, Grass & Heckel**, "Embracing errors is more efficient than avoiding them through constrained coding for DNA data storage". [arXiv:2308.05952](https://arxiv.org/abs/2308.05952); IEEE TMBMC | Argues that randomisation plus more ECC beats constrained coding unless problematic sequences have very high error rates. | **Directly challenges our premise.** See the overlap warning, item 4. |
| **Gungnir (2026)**, *Nat. Commun.* [doi:10.1038/s41467-026-71485-x](https://doi.org/10.1038/s41467-026-71485-x) | Hash-signature, guess-until-correct decoding; recovery from one copy with 20% errors; considers GC, homopolymer, and error-prone motifs; aimed at nanopore. | Recent state of the art on the decoding side. It shows that compute-heavy *decoding* can replace constraint-heavy *encoding*. |

## 3. Secondary structure / MFE

| Work | What it does | Overlap |
|---|---|---|
| **Cao et al. (2021)**, "Minimum free energy coding for DNA storage", *IEEE Trans. NanoBiosci.* 20(2):212–222. [doi:10.1109/TNB.2021.3056351](https://doi.org/10.1109/TNB.2021.3056351) | Adds an MFE constraint when constructing DNA coding sets. | MFE as an encoding objective already exists. Our contribution is limited to MFE *inside a per-oligo search, including primers and both strands*. |
| **Nguyen et al.**, secondary structure avoidance (SSA) codes. [arXiv:2302.13714](https://arxiv.org/abs/2302.13714); improved constructions [arXiv:2304.11403](https://arxiv.org/abs/2304.11403) | Combinatorial codes that forbid reverse-complement stems of length ≥ m (rate 1.3031 bits/nt for m = 3). | Gives a guarantee where we use a heuristic. Our MFE scoring is thermodynamic, which is more realistic, but it comes with no guarantee. |

## 4. Nanopore-aware and signal-level coding

| Work | What it does | Overlap |
|---|---|---|
| **Whritenour, Civelek & Farnoud (2025)**, "Constrained coding for error mitigation in nanopore-based DNA data storage", *Sci. Rep.* [doi:10.1038/s41598-025-08531-z](https://doi.org/10.1038/s41598-025-08531-z) *(read via PMC summary)* | Adjacent k-mer current difference > τ (R9.4.1, 6-mer, ONT tables), RLL ≤ 3, state-splitting encoder; DeepSimulator signal; Guppy/Dorado/Viterbi. About 6× lower edit distance at τ = 4 pA with a rate cost of 0.7 bits/base. | **Closest prior work to the confusability objective.** See the overlap warning, item 1. |
| **Vidal, Wijekoon & Viterbo (2024)**, "Concatenated nanopore DNA codes", *IEEE Trans. NanoBiosci.* 23(2):310–318. [doi:10.1109/TNB.2024.3350001](https://doi.org/10.1109/TNB.2024.3350001) *(title and venue from search; not read)* | Codes designed for a nanopore channel model. | Nanopore-specific coding. To be read. |
| **"On coding for an abstracted nanopore channel for DNA storage"**. [arXiv:2102.01839](https://arxiv.org/abs/2102.01839) *(not read)* | Information-theoretic nanopore channel abstraction (k-mer ISI). | Theory for channel models in which k-mer confusability matters. |
| **Chandak et al. (2020)**, basecaller–decoder integration with convolutional codes, ICASSP 2020. [bioRxiv 10.1101/2019.12.20.871939](https://doi.org/10.1101/2019.12.20.871939) | Viterbi decoding on basecaller soft output; about 3× lower reading cost. | Decodes at signal level instead of shaping sequences. Complementary to our approach. |
| **Doroschak et al. (2020)**, Porcupine. [bioRxiv 10.1101/2020.03.06.981514](https://doi.org/10.1101/2020.03.06.981514); published in *Nat. Commun.* as "Rapid and robust assembly and decoding of molecular tags with DNA-based nanopore signatures" | 96 "nanopore-orthogonal" molbits designed to be separable in **raw signal space**, classified without basecalling. | Designs sequences for distinguishability in signal space, which is the same motivation as our confusability score (applied to tags, not payload). |
| **Composite Hedges Nanopores (2024)**, *Nat. Commun.* [doi:10.1038/s41467-024-53455-3](https://doi.org/10.1038/s41467-024-53455-3) | HEDGES-derived codec for nanopore readout; tolerates up to 15.9% indels and 7.8% substitutions. | Nanopore-targeted codec, but on the ECC side rather than sequence shaping. |

## 5. JPEG DNA and image-specific codecs

| Work | What it does | Overlap |
|---|---|---|
| **JPEG DNA, ISO/IEC 25508-1** (DIS stage, 2026). [jpeg.org/jpegdna](https://jpeg.org/jpegdna/); [ISO 90579](https://www.iso.org/standard/90579.html) | Standardised image coding into constrained quaternary representations. | Different philosophy: source-aware (works on pixels/coefficients). Our spec deliberately encodes **file bytes** and treats media as opaque. Note this contrast in the paper. |
| **Dimopoulou, Antonini et al. (2019)**, "A biologically constrained encoding solution for long-term storage of images onto synthetic DNA", EUSIPCO 2019. [arXiv:1904.03024](https://arxiv.org/abs/1904.03024) | DWT + quantisation + a constrained quaternary code. | Image-specific constrained coding. |
| **"A JPEG-based image coding solution for data storage on DNA"**. [arXiv:2103.09616](https://arxiv.org/abs/2103.09616); **"A constrained Shannon-Fano entropy coder for image storage in synthetic DNA"**. [arXiv:2203.09988](https://arxiv.org/abs/2203.09988) | Precursors of JPEG DNA. | Same as above. |

## 6. Epigenetic and non-sequence bit storage

| Work | What it does | Overlap |
|---|---|---|
| **Zhang et al. (2024)**, "Parallel molecular data storage by printing epigenetic bits on DNA", *Nature* 634:824–832. [doi:10.1038/s41586-024-08040-5](https://doi.org/10.1038/s41586-024-08040-5) | Writes bits as methylation on premade templates by self-assembly-guided enzymatic methylation; nanopore readout of the modifications. | Orthogonal: a different write channel. Relevant only because nanopore methylation calling is another signal-level readout. No codec overlap. |
| **Tabatabaei et al. (2020)**, "DNA punch cards", *Nat. Commun.* 11:1742. [doi:10.1038/s41467-020-15588-z](https://doi.org/10.1038/s41467-020-15588-z) | Data written as nicks on native dsDNA. | Orthogonal. |

---

## Gaps in this scan (to do before writing the paper)
- No systematic search for "redundant / pilot / flexible base" schemes from 2023–2026. Several Chinese-group codecs (e.g. DNA-MGC+ [arXiv:2603.14527](https://arxiv.org/abs/2603.14527)) were not examined.
- I did not read the YYC full text to confirm whether it screens on free energy.
- I did not look at Illumina-specific error-motif literature (e.g. GGC motifs), which would bear on a "readout-aware" claim beyond nanopore.
