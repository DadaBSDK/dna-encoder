# Read-channel calibration: real R10.4.1, squigulator and Badread

Phase 3 calibration, 2026-10-03. Reproduce with `scripts/calibrate_channel.py`; raw calls,
alignments, summaries and reports are in `results/phase3/calibration/`.

## Dataset and comparison

The real signal is the public Chen et al. 2025 DNA-storage dataset (MinION FLO-MIN114,
SQK-RBK114.24, 5 kHz; Zenodo record [10.5281/zenodo.16883332](https://doi.org/10.5281/zenodo.16883332),
CC-BY-4.0; also SRA PRJNA1235219). It contains 33 plasmids with encoded inserts, not a pool
of 200-nt storage oligos. Plasmids are 6–43 kb and median read length is about 355 nt. The
comparison uses identical reference intervals and strands, Dorado hac v6.0.0, MAPQ ≥ 20,
105 real reads (1.1 Mb) and mean scored depth about 7.8×. The encoded-region run-length
analysis is separated from plasmid-backbone sequence. MAPQ ≥ 0 sensitivity analysis uses
335 reads and gives 1.55% hac and 1.05% sup per-base error, close to the MAPQ ≥ 20 result.

Badread v0.4.1 uses its nanopore2023 error model and identity parameters
`98.48,99.9,1.2`, chosen to match the mean and standard deviation of real hac per-read
identity. It is an independent read-level channel, not a signal simulator. It does not
model the read-to-read systematic errors seen in real reads.

## Results

| Measure | Real hac | Real sup | Squigulator → Dorado hac | Badread, identity-matched |
|---|---:|---:|---:|---:|
| Per-base error | **1.51%** | 1.01% | 9.12% | 1.55% |
| Substitution / deletion / insertion | 0.67 / 0.57 / 0.27% | — | 4.7 / 3.0 / 1.4% | 0.67 / 0.58 / 0.30% |
| Encoded-region error by run length 1 / 2 / 3 / 4 / 5 / 6+ | 1.56 / 1.31 / 1.40 / 1.75 / 2.76 / 3.36% | — | 9.1 / 8.3 / 10.4 / 12.3 / 16.6 / 22.1% | 1.5 / 1.2 / 1.5 / 1.7 / 1.7 / 3.7% |
| Systematic-site fraction, run 1–3 / 4 / 5 / 6+ | ~0.01 / 0.06 / 0.35 / 0.66% | ≤0.11% | 2.9–4.3 / 6.0 / 11.7 / 18.8% | 0 |

Systematic sites are reference positions where at least half of covering reads contain an
error. The real encoded sequence retains a homopolymer gradient, with a detectable systematic
component mainly at runs ≥ 5, but its absolute rate is small. Squigulator produces about six
times the raw error and about 300 times the short-run systematic error. This is consistent
with its documented crude 5 kHz R10 model being mismatched to Dorado's training signal.

## Evaluation implications

- Use identity-matched Badread as the practical read-level channel for broad end-to-end
  comparisons, and retain uniform IDS as a controlled channel.
- Keep squigulator → Dorado as a **stress test**, clearly labelled miscalibrated. Its
  1.5–2.4% per-base consensus floor across 30–40 reads leaves only 4–7% of oligos within
  the current two-edit run-length repair radius; no arm recovers a file in that channel.
- Badread is not a replacement for real reads: its independent per-read errors omit the
  systematic sites observed in the real data. The dataset is plasmid-fragment data, not a
  storage-oligo pool, so transfer to 200-nt pools is approximate.
- Any optional systematic homopolymer stress term must be reported separately from the
  identity-matched Badread result. Do not describe the Badread channel alone as modelling
  real consensus floors.

The later binding DESIGN.md §0c.8 adopts Badread as the main read-level evaluation,
uniform IDS as the controlled comparison, and squigulator as a stress channel. The
matched-budget benchmark follows that decision; the limitations above still apply.
