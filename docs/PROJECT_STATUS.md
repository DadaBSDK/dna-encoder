# Project takeover and status — 2026-10-08

## Scope and review

Read the project README, binding DESIGN revisions, RELATED_WORK and all project-owned
methodology/handoff Markdown files. Reviewed the encoder, decoder, ECC/container layout,
channel and CLI implementations, existing tests and experiment outputs. The supplied
folder has no Git metadata. Existing research assets and the already-running 100 KB
clustering experiment were preserved.

The current research plan already selects calibrated Badread as the main read-level
channel, uniform IDS as the controlled comparison and squigulator/Dorado as a
miscalibrated stress channel. Older notes presenting that choice as unresolved were
updated. No literature claims were newly validated in this software review.

## Completed implementation

| Finding | Change |
|---|---|
| DESIGN promised a benchmark command, but none existed | Added `dnastore benchmark`, exact physical nucleotide matching, independent channel seeds, both decoder strengths, Wilson CIs, empirical integer depth bisection, CSVs, plots and provenance |
| Default primer lookup depended on the current directory | Bundled default primers in the wheel; verified byte-exact recovery from outside the checkout |
| CLI silently exempted baselines from the strict hard-constraint policy | All codecs now require `--allow-violations` to emit violating pools; diagnostic manifests still record violations |
| FASTQ parsing assumed four lines and ignored malformed records | Added wrapped FASTQ, multiline FASTA and gzip handling with explicit length/format checks |
| Reads containing `N` could reach base lookup code and crash | Ambiguous reads are discarded and counted; lowercase ACGT reads are normalized |
| Invalid RS, seed, header-copy, channel and codec settings failed late or behaved nonsensically | Added range/type checks and clear CLI errors |
| Header parsing accepted non-object codec parameters and some invalid geometry | Added header-parameter/geometry checks and structured decode failures |
| Repair mode attempted data repairs on valid header-only clusters | Skip groups only when every candidate is a CRC-verified header; preserve uncertain and mixed groups |
| Decoder strengths were inaccessible through the CLI | Added `--strength erasure|repair` and `--grouping cluster|index`; CLI defaults to repair |
| No automated project workflow | Added a GitHub Actions test/build workflow, using pinned model data; not executed on GitHub because this folder has no remote |
| Stale channel guidance and missing operating instructions | Updated README and method notes; added `docs/BENCHMARK.md` and smoke/study configurations |

The on-disk format remains version 3. Existing valid pools remain compatible. Deliberate
behavior changes are stricter input validation and the CLI baseline constraint policy.
`decode_reads` retains its existing API default decoder strength.

## Verification

Initial baseline: **205 passed, 1 skipped**. The expanded full suite passed with
**240 passed, 1 skipped**. After the final header-repair optimization, the affected
CLI/benchmark/repair/noisy-decoder suites passed **44 tests**, including its new regression.
Both header-only and mixed-header/data-cluster regressions also passed separately.
The pre-existing skip covers the alphabet-C steering semantic test whose
behavior is covered by a separate C-specific test.

- Dependency consistency: `pip check` passed.
- Package wheel built successfully using the declared isolated build environment.
- Extracted-wheel import, bundled primer loading and byte-exact roundtrip passed with
  `/tmp` as the working directory, independently of the checkout import path.
- All modified Python modules compile.
- Six-arm IDS smoke run: `results/phase4/smoke_20261008/`, **152 decoder trials**,
  exactly **62 × 200 = 12,400 nt per arm**. All six noiseless checks succeeded, and all
  six arms recovered both channel seeds at depth 8 with both strengths. Two-seed depth
  boundaries are descriptive pipeline checks, not scientific evidence.
- Calibrated Badread smoke: `results/phase4/badread_smoke_20261008/`, **8 decoder trials**,
  whiten-RS and Goldman, exactly **49 × 200 = 9,800 nt per arm**, both strengths and both
  channel seeds recovered at requested depth 8.
- PNG recovery plots inspected; SVG equivalents and source CSVs are saved alongside them.
- Multi-seed validation: `configs/benchmark_validation.yaml` / `results/phase4/validation_optimized_20261008/`.
  Two 800-byte payload seeds, six arms, 20 independent IDS channel seeds at depth 8,
  both strengths: **480 trials completed**, all pools passed noiseless recovery. Each arm
  used 62 oligos (12,400 nt). Whiten-RS, Goldman, steering-A-P8, steering-C-P8 and fountain
  recovered 20/20 seeds at both strengths for both payloads. Steering-B-P6 recovered
  17/20 (erasure) and 20/20 (repair) on payload 0, and 19/20 at both strengths on payload 1.
  These are small synthetic-payload validation results, not a full corpus comparison.
  Each payload's Wilson intervals are in `summary.csv`. Single-depth figures were later
  redrawn with separate codec positions for readability; `plot_provenance.json` records
  the unchanged trial-data hash and the plotting-module hash.
- Profiled the same saved read realization before/after the header-only optimization:
  17 → 5 data-repair attempts, ~32,000 → ~8,000 generated candidates, 9.32 → 3.53 s
  under cProfile with concurrent background load. Both runs recovered identical bytes.
  Across 120 shared channel/strength trials, pre/post-optimization recovery and
  continuous error metrics were identical (runtimes excluded).
  This is a single-case performance check, not a general throughput claim. The first
  multi-seed run (`validation_20261008`) was explicitly interrupted and preserved when
  this issue was found; the optimized run uses a new directory and source snapshot.

## Research work that remains

The software now has a working experiment runner. The research study itself is not
finished, and this review does not claim a winning codec or publication-ready evidence.

The existing **100 KB clustering scale run finished during this review**: all 12
arm/channel rows are present, and the process has exited. All six arms recovered under
Badread at both strengths; none recovered under IDS9 at depth 10. All rows have zero
cluster merges. Badread split rates are zero; IDS9 split rates range from 2.04% (fountain)
to 6.41% (steering B). Minimum distinct-body distances range from 0.31875L (steering B)
to 0.38125L (whitened quaternary); no pair lies within the 0.30L cluster radius.
Here IDS9 uses 4% substitutions + 2.5% insertions + 2.5% deletions, as specified in
`cluster_scale.py`. This scale run uses one payload and unmatched nucleotide budgets;
it checks clustering behavior, not whether constraints outperform matched ECC.
Source: `results/phase3/cluster_scale_100KB/cluster_scale.csv`.

1. **Run the larger Phase 4 design.** `configs/paper.yaml` is a runnable core study
   preset (three random 4 KB payloads, seven arms, three channels, at least 20 seeds,
   depth search). Its pool matcher was preflighted successfully at 245 oligos per arm
   on an incompressible 4 KB payload. It has not been run to completion in this review.
2. **Add the publication corpus and parameter sweeps.** Real, licensed files with URLs
   and hashes; steering periods and weight-grid/ablations; the 100 KB and larger timing
   cases; confusability delta sweeps and adjacent-level contrasts using the existing
   characterization scripts. Use pairwise matched runs when a steering setting should
   determine its own budget, as described in BENCHMARK.md.
3. **Keep model limitations explicit.** Badread omits systematic inter-read errors;
   calibration uses plasmid fragments. Accessibility-to-yield remains an assumed
   sensitivity parameter, with paired on/off results required. The main benchmark
   deliberately rejects circular and assumed-accessibility channel settings.
4. **Publication preparation.** Complete the literature verification gaps already
   listed in RELATED_WORK.md, analyze the full study, and write conclusions only from
   those results. Historical Phase 2 density tables predate CRC-32 and should not be
   mixed directly with version-3 benchmark measurements.

Long experiments save completed depth probes incrementally and mark failed/interrupted
runs in provenance. Automatic resume is not implemented; keep old output directories and
use a new one for each run. Source archives are included in runs launched after the
initial smoke checks, so later work can recover the implementation without Git.


## Streaming media extension — 2026-10-08

Added `stream-encode` and `stream-decode` for arbitrary file types, including images,
animated GIFs and videos. This separate version-1 DNA text archive processes bounded
chunks without a configured total nucleotide or file-size cap. It preserves bytes,
records the original filename, and verifies SHA-256 and length before publishing
recovered output. Existing outputs are protected. The existing version-3 physical
oligo format and its error correction remain unchanged.

This digital archive mode provides integrity detection, not sequencing-error correction
or synthesis screening. DNA payload text is four times the input size, plus framing.
Actual media demonstrations and hashes are in `results/media_stream_examples/results.json`;
`scripts/media_examples.py` reproduces them. See README for usage and format boundaries.

Verification after this extension: full suite **253 passed, 1 skipped** (185.35 s);
`pip check` passed. Valid PNG (8,975 bytes), animated GIF (17,020 bytes), and MP4
(410,964 bytes) each recovered exactly and decoded successfully with FFmpeg. The MP4
produced 1,643,856 DNA payload nucleotides. These are digital round-trip demonstrations,
not noisy-sequencing recovery measurements.
