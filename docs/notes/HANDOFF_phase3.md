# Phase 3 handoff (written 2026-10-03; read this first when resuming)

Historical handoff. Current implementation and remaining work are tracked in
[`../PROJECT_STATUS.md`](../PROJECT_STATUS.md). The later binding DESIGN.md §0c.8
resolves the channel decision discussed below.

## User instructions in force (Phase 3 revision, from the user)
1. Decoder: CRC-guided repair (homopolymer run-length ±1, ≤ 2 edits, accept on length + CRC), inner CRC-16 → CRC-32, outer RS errors-and-erasures. One decoder for all arms. Report **two decoder strengths**: `erasure` and `repair`.
2. Keep the uniform IDS channel as a controlled second channel for file recovery and min-depth, alongside squigulator → Dorado.
3. Calibrate the simulator against real public R10.4.1 storage data, and evaluate Badread as an independent channel. If real data shows a much weaker homopolymer effect, say so and rescale.
4. Clustering scale check at Phase 4 sizes (thousands of oligos): min distinct-oligo distance, merge rate, and how often sub-grouping splits groups.
5. Add the sweep_v2 results (with alphabet C) to DESIGN.md and send the user the updated characterisation table.

Then continue Phase 3. The standing rule applies: STOP after each phase, summarise, and wait. Never fabricate.

## Status per item
| # | Status |
|---|---|
| 1 | **DONE (code + tests).** `dnastore/repair.py`; `ecc.py`: CRC-32, `FORMAT_VERSION = 3`, errors-and-erasures RS; `decoder.py`: strengths `erasure` and `repair`, CRC-guided repair and reports. Current full suite: **205 passed, 1 skipped** (2026-10-08). |
| 2 | **DONE.** `scripts/nanopore_e2e.py --channel ids` evaluates uniform IDS through grouping, depth subsampling and both decoder strengths. `--ids-rate` defaults to 9% total events split equally; an 800-byte smoke run completed for whitened 2-bit and Goldman (no file recovery at the smoke depth). |
| 3 | **DONE.** Measurement and documentation: `docs/notes/channel_calibration.md`; `docs/notes/nanopore_pipeline.md` now records the systematic squigulator floor. Evaluation plan: Badread (main read-level), IDS (controlled), squigulator (stress). |
| 4 | **IN PROGRESS.** The 8 KB smoke run finished; results in `results/phase3/cluster_scale_smoke8KB/cluster_scale.csv`. The requested 100 KB, six-arm job is running under `results/phase3/cluster_scale_100KB/`:
  - **Min distance already shrinks at 300–400 oligos:** 0.356·L (whiten), **0.331·L (goldman)**, compared with 0.36 at ~245 oligos. Zero pairs ≤ 0.30·L. 19 goldman pairs ≤ 0.36·L.
  - **Clustering:** merge rate 0 in both channels. Split rate 0 on Badread, ~2% of oligos on IDS-9%. `groups_multi_index` = 0 everywhere.
  - **Decoding at 10×:** Badread decodes at both strengths (repair fixed 5 oligos for whiten). IDS-9% fails at both strengths; repair raises valid oligos from 170 to 214 for whiten, 262 to 265 for goldman.
  - **Current run:** `.venv/bin/python scripts/cluster_scale.py --bytes 100000 --out results/phase3/cluster_scale_100KB`; inspect `cluster_scale.csv` when it exits. |
| 5 | **DONE.** Added the alphabet C sweep and full characterization table to `docs/DESIGN.md` §0c.7, sourced from `results/phase2/sweep_v2/summary.csv`. |

## Key findings so far
**Decoder issues found and fixed (Phase 3, part 1):**
- The primer search window was 12 nt, which rejected all adapter-flanked reads; it is now 80 nt.
- Exact-index grouping fails at ~10% read error (26% of reads get the correct index). The decoder now uses `kmer_cluster` (k = 10, read-to-read bound 0.30·L) plus sub-groups by per-read index.
- Measured distances (~245 oligos): same-oligo Dorado read pairs median 0.156·L, p99 0.287, max 0.33. Distinct oligos min 0.36·L.

**squigulator → Dorado hac has a systematic consensus floor:** 1.5–2.4%/nt after consensus over 30–40 reads. It holds for both the dna-r10-prom and dna-r10-min profiles, hac and sup. By arm: whiten 2-bit 3.9 edits per 160 nt, fountain 2.5, goldman 1.9. No arm decodes a file at any depth or strength.
- **Repair ceiling on squigulator** (`results/phase3/e2e_diag/repairability.py`): only 4–7% of consensus bodies are within 2 run-length edits of the truth. 35–52% of the consensus errors are substitutions.

**Calibration result** (`scripts/calibrate_channel.py`; `results/phase3/calibration/`):
- Real data: Chen et al. 2025, Nat. Commun., doi:10.1038/s41467-025-65004-7. Zenodo 10.5281/zenodo.16883332 (CC-BY-4.0; zip sha256 `8006d9c9…63c5`), stored in `data/real_r10/chen2025/`. FASTQ also in SRA PRJNA1235219.
- The run: MinION FLO-MIN114, SQK-RBK114.24, 5 kHz, 33 plasmids (6–43 kb) with encoded inserts. **Caveat:** these are plasmid fragments (median read 355 nt), not a 200-nt oligo pool.
- Also found: a Helixworks R10.4.1 composite-motif dataset (Zenodo record 15839178, FAST5). Not used.
- Comparison: identical reference intervals and strands, Dorado hac v6.0.0, MAPQ ≥ 20 (105 reads, 1.1 Mb), depth ≈ 7.8× at scored positions:

| | real hac | real sup | squigulator hac | Badread (identity matched) |
|---|---|---|---|---|
| per-base error | **1.51%** | 1.01% | 9.12% | 1.55% |
| sub / del / ins | 0.67 / 0.57 / 0.27% | — | 4.7 / 3.0 / 1.4% | 0.67 / 0.58 / 0.30% |
| error, encoded, run 1/2/3/4/5/6+ | 1.56/1.31/1.40/1.75/2.76/3.36% | — | 9.1/8.3/10.4/12.3/16.6/22.1% | 1.5/1.2/1.5/1.7/1.7/3.7% |
| systematic sites (≥ 50% of reads err), encoded, run 1–3 / 4 / 5 / 6+ | **~0.01% / 0.06% / 0.35% / 0.66%** | ≤ 0.11% | 2.9–4.3% / 6.0% / 11.7% / 18.8% | 0 |

- MAPQ ≥ 0 sensitivity (335 reads): 1.55% (hac), 1.05% (sup). Robust.
- Real per-read identity: mean 98.48, sd 1.17. Badread default set to `98.48,99.9,1.2`.
- **Conclusion:** squigulator's floor is a simulator artefact, about 300× the real systematic rate at short runs. Real data keeps a homopolymer gradient (systematic component only at runs ≥ 5), at tiny absolute levels. Badread matches the real per-base and per-run-length rates but has **no** systematic component.
- Badread smoke run (800 B, ~40×): whiten-2bit and goldman decode at 10×. Fountain fails even noiselessly at this size (LT with k = 26 at 15% overhead is a small-k issue, not a bug); use ≥ 4 KB.

## Decision to put to the user when reporting
squigulator → Dorado was the "required" main channel, but it is miscalibrated: 6× the raw error and ~300× the systematic error. Proposal:
- **Main:** Badread (calibrated), plus an optional calibrated systematic homopolymer term for runs ≥ 5.
- **Controlled:** uniform IDS.
- **Stress / appendix:** squigulator, clearly labelled.

The user decides.

## Other changes this session
- `tests/test_consensus.py`: flanked-read test. `tests/test_fountain_codec.py`: tiny-k LT may fail but must fail loudly.
- Installed in `.venv`: mappy 2.31, Badread v0.4.1 (git tag).
- `scripts/nanopore_e2e.py`: `--channel`, `--strengths`, per-oligo error including the present-only metric, `oligos.fa` saved per run.
- Memory file updated: `~/.claude/projects/-home-devashish-wakde-Desktop-DNA-Encoder/memory/dna-codec-project.md`.
