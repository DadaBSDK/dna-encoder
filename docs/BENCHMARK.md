# Reproducible matched-budget benchmark

Implemented 2026-10-08 in `dnastore/benchmark.py`; exposed as `dnastore benchmark`.
This completes the experiment runner. It does not by itself complete the publication's
corpus, period, weight-grid, or sensitivity experiments in `DESIGN.md`.

## Run

From the project root, after installation:

```sh
.venv/bin/dnastore benchmark configs/benchmark_smoke.yaml --out results/phase4/my-smoke
.venv/bin/dnastore benchmark configs/paper.yaml --out results/phase4/my-study
```

The output directory must be new or empty. Existing results are never overwritten.
The smoke preset uses six arms, 800 random bytes, two independent channel seeds and
uniform 3% IDS. It exercises actual default steering weights, normalization and
accessibility rescoring. **Smoke outcomes cannot establish a scientific advantage.**
The study preset uses three 4 KB payload seeds, seven arms, calibrated Badread plus
3% and 9% IDS, and 20 independent channel seeds at every probed depth. Both decoder
strengths are always evaluated on the same reads. This can take hours on a laptop.
The runner saves completed depth probes incrementally; automatic resume is not implemented.

Badread is optional and required only for `kind: badread`. The existing integration
expects it beside the running Python executable (`.venv/bin/badread` in this setup).
Install the pinned upstream tag:

```sh
.venv/bin/python -m pip install 'git+https://github.com/rrwick/Badread.git@v0.4.1'
```

Steering's default confusability objective requires the model tables. If they are missing,
run `.venv/bin/python scripts/fetch_kmer_models.py`. Source checksums are verified by the
model loader. Badread calibration and its limits are in `notes/channel_calibration.md`.
Raw-signal stress tests remain in `scripts/nanopore_e2e.py --channel squigulator`.

## Exact nucleotide matching

For each payload, every arm receives the identical stored byte stream, oligo length,
primer pair, seed-field width and number of copies per header fragment. The matcher
searches integer RS parity ratios from the configured minimum (150 per thousand by
default) through `max_parity_permille` (2000 by default). It chooses the smallest pool
size attainable by **every** RS arm; at that size each arm uses the lowest matching
ratio. Block size follows the existing maximum-valid GF(256) rule, unless explicitly
fixed in the encode configuration.

All physical oligos count: header fragments and their copies, indices, seeds, CRC,
steering, parity, primers and padding. Header fragment counts differ by codec metadata.
The fountain arm receives the exact remaining number of droplets after its own header
cost. Its configured overhead is a lower bound. No unused filler oligos are added.
The actual pool sizes are asserted after encoding. Impossible budgets fail explicitly.

This global matching gives all arms a common budget. For a pairwise P-sweep, make one
configuration per steering setting containing that arm and whiten-RS; otherwise a
low-rate arm can set a much larger budget for every other arm. Record the matched ratios
in `pools/*/*/manifest.json`. Very small payloads magnify header overhead; do not interpret
the smoke preset's density as representative of large files.

## Configuration

Paths in `encode_config` and payload `path` are relative to the benchmark YAML.
The encode configuration resolves its own primer path relative to its own directory.
Without `encode_config`, the package's default settings and bundled primers are used.

| Field | Meaning |
|---|---|
| `mode` | `study` (default, at least 20 channel seeds) or `smoke` |
| `seed` | Root nonnegative RNG seed |
| `encode_config`, `encode` | Base encoder YAML and optional nested overrides |
| `payloads` | Named entries, each with either positive `bytes` (synthetic random data) or `path` |
| `payload_seeds` | Unique nonnegative seeds; default `[0]`. Files stay identical, encoder randomization changes. |
| `channel_seeds` | Independent full channel realizations per depth, default 20 |
| `arms` | Unique named entries with `codec: {name, params}` |
| `channels` | Unique named entries with `kind: ids` or `badread`, and `params` |
| `depths` | Unique positive integer mean depths, default `[3, 5, 10]` |
| `max_depth` | Upper bound for depth search, at least the largest requested depth |
| `minimum_depth` | Enable empirical integer bisection (default true) |
| `max_parity_permille` | Upper bound for budget matching, default 2000 |
| `plots` | Export PNG and SVG recovery plots (default true) |

IDS parameters are `ChannelParams` fields, with coverage and seed controlled by the
runner. Badread accepts `identity`, coverage parameters and dropout; the model itself
sets error rates and strand mixture. Actual observed depth is saved because Badread's
base-count target only approximately produces the requested number of reads.
Circular confusability-weighted and assumed-accessibility coverage models are excluded
from this main runner; use the existing sensitivity modules with their required labels
and paired accessibility-on/off runs.

## Statistics and reproducibility

Each depth gets independent channel realizations, rather than repeated subsamples of a
single simulated library. Random streams derive from `numpy.random.SeedSequence` using
stable, SHA-256-derived labels. Reordering arms does not change seeds; both strengths see
identical reads, and matched arms share seed labels. Channel seeds are independent
within each payload/depth; seeds do not remove the need for multiple source files.

`summary.csv` reports recovery fractions and Wilson 95% confidence intervals **separately
for each payload seed**, avoiding pseudoreplication across shared payloads. Continuous
metrics include recovered/attempted-body edit distance per nt, missing oligos and the
fraction not recovered with valid CRC before outer decoding. Missing bodies score 1.0
in `body_error_mean`; `body_error_present_mean` in the trial table excludes them. These
are decoder-output metrics and can include CRC-guided repair; they are not a measurement
of a separate, untouched consensus stage. The unweighted mean gives each oligo one vote.

Depth search probes the upper bound and bisects for the first observed all-seed success.
Its result is an **empirical bisection boundary**, not a guaranteed physical minimum:
independently sampled recovery outcomes need not be monotone. `minimum_depth.csv` also
records the lowest observed all-success depth, right censoring at the upper bound, and
nonmonotonicity among the probed depths. Every probe's Wilson interval remains in the
summary. Even 20/20 successes has a lower 95% bound of about 0.839, not 100% certainty.

A noiseless decode is recorded for every pool. A rank-deficient fountain pool is a
measured failure, never silently retried or excluded. Every successful trial requires
both decoder integrity checks and exact byte equality with the original input.

## Output files

- `config.yaml`, `encode_config.json`: requested settings and resolved base settings.
- `provenance.json`: completion/failure state, source hashes, package versions, primer
  and payload hashes, root seed, k-mer provenance, timestamps and Git SHA when available.
  This supplied project folder has no Git metadata; hashes provide the source identity.
- `source.zip`: a snapshot of the Python implementation, bundled primers and package
  metadata, preserving the implementation even without Git. Added after the initial
  smoke runs; subsequent runs include it automatically.
- `encode.csv`: actual nucleotide costs, parity, density, runtime, constraint violations
  and noiseless recovery for every arm/payload.
- `pools/<payload_seed>/<arm>/`: FASTA and detailed configuration/screening manifest.
- `trials.csv`, `summary.csv`: individual channel trials and conditional recovery CIs.
- `minimum_depth.csv`: empirical depth boundaries and censoring (when enabled).
- `recovery_*.png`, `recovery_*.svg`: recovery curves with Wilson intervals.
- `channels/`: Badread references, commands and generated reads when that channel is used.

The manifest records best-effort constraint violations for all benchmark arms. In the
ordinary `encode` CLI, any hard-constraint violation requires `--allow-violations`,
including for baseline codecs, as specified in DESIGN.md D4.
