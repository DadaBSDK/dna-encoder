# DNA storage codec research project

Encode file bytes into fixed-length DNA oligos and recover them from sequencing reads.
Shared indexing, CRC-32 and Reed–Solomon error correction support quaternary, Goldman
rotating ternary, steering A/B/C and LT fountain payload mappings.

The research question is whether sequence constraints improve recovery enough to justify
the nucleotides they consume. The binding design and prior findings are in
[`docs/DESIGN.md`](docs/DESIGN.md); scope and overlap with prior work are in
[`docs/RELATED_WORK.md`](docs/RELATED_WORK.md). This is research software, not a
validated synthesis or clinical workflow.

## Setup

Python 3.11 or newer:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/dnastore --help
```

Core encoding, decoding and IDS simulation need no raw-signal tools. Default primers
are bundled with the package, so the CLI also works outside this folder. Explicit
configuration files resolve primer paths relative to their own directory.

Default steering uses ONT confusability tables. Fetch them if absent:

```sh
.venv/bin/python scripts/fetch_kmer_models.py
```

Optional Badread and raw-signal dependencies are documented in
[`docs/BENCHMARK.md`](docs/BENCHMARK.md) and
[`docs/notes/nanopore_pipeline.md`](docs/notes/nanopore_pipeline.md).

## Encode and recover a file

```sh
.venv/bin/dnastore encode README.md --codec goldman --allow-violations \
  --out-dir results/example
.venv/bin/dnastore decode results/example/README.md.goldman.fasta \
  --out results/example/recovered --strength repair
```

The decoder chooses the recovered extension from the encoded file type (`.txt` for
this README). Successful recovery requires the original length and file CRC-32 to match.
The format is version 3; keep pools and decoders on compatible format versions.

Encoding writes a manifest with configuration, density and constraint violations.
Without `--allow-violations`, **all codecs** refuse to write a pool with hard-constraint
violations. The option above is explicit because baseline mappings can violate sequence
constraints. It does not bypass integrity checks.

Input reads may be multiline FASTA, wrapped FASTQ, or either format compressed as `.gz`.
Malformed records produce an error. Reads containing `N` are discarded and counted in
`ambiguous_reads`. `--strength erasure` selects CRC-discard/RS-erasure decoding;
`--strength repair` (CLI default) also enables bounded CRC-guided homopolymer repair
and RS errors-and-erasures. API defaults remain unchanged.

## Images, videos, GIFs and uncapped DNA streams

Both modes accept arbitrary file bytes, including images, videos and animated GIFs.
They preserve the encoded bytes, so resolution, frames, sound and metadata are unchanged
when recovery succeeds. The oligo container already recognizes PNG, JPEG, GIF, WebP
and MP4; unknown formats remain supported as bytes and recover with `.bin` in that mode.

For large digital media, use the **streaming mode**:

```sh
.venv/bin/dnastore stream-encode video.mp4 --out video.mp4.dna
.venv/bin/dnastore stream-decode video.mp4.dna --out recovered-video.mp4
```

The same commands work for `photo.png`, `animation.gif`, or any other file type.
There is **no configured maximum file size or total nucleotide count**. Input is read
in 64 KiB chunks, so memory does not grow with the entire file. Disk space and runtime
remain practical limits. Each byte maps to four nucleotides (`00=A`, `01=C`, `10=G`,
`11=T`); DNA text takes approximately four times the input space plus line/framing
metadata. 1 GiB therefore needs about 4 GiB of DNA text. This is encoding, not compression.

A `.dna` archive contains a versioned metadata header, DNA lines, and a SHA-256/length
footer. The original filename is recorded; `--out` explicitly chooses the recovered
filename including extension. Lines wrap at 1,024 nucleotides for bounded processing;
this is a formatting choice, **not a total sequence limit or physical oligo length**.
The header and footer are text metadata, so the archive is not a synthesis-ready FASTA.

Streaming mode detects corruption but does **not** correct mutations, missing reads,
or reordering, and does not screen synthesis constraints. Use the existing `encode` /
`decode` oligo workflow for error-corrected sequencing experiments. Its format-specific
index and header bounds still apply; changing `oligo_len` does not remove those bounds.

Both streaming commands refuse to overwrite existing output. The decoder publishes
recovered output only after validating the archive length and SHA-256. Failed operations
remove temporary files. SHA-256 detects accidental corruption; this format is not an
authenticated or encrypted archive.

Try streaming mode without installing anything at <https://dadabsdk.github.io/dna-encoder/>. The page (`site/`) is a browser port of the same
format: archives are byte-identical to `stream-encode`, decoding performs the same checks, and files
never leave the browser. `tests/test_site_stream.py` checks interoperability with Node.js.

Reproduce valid PNG, animated GIF and MP4 demonstrations (FFmpeg needed only to create
and validate the fixtures; the output directory must not already exist):

```sh
.venv/bin/python scripts/media_examples.py
```

## Simulate a read channel

```sh
.venv/bin/dnastore simulate results/example/README.md.goldman.fasta \
  --out results/example/reads.fasta --seed 42
.venv/bin/dnastore decode results/example/reads.fasta \
  --out results/example/from_reads
```

The default is error-free reads at mean depth 10 with random strand orientation and
multinomial coverage. Use `--channel path/to/channel.yaml` to set `p_sub`, `p_ins`,
`p_del`, dropout and coverage parameters. All probabilities are validated.

## Matched-budget benchmark

```sh
.venv/bin/dnastore benchmark configs/benchmark_smoke.yaml \
  --out results/phase4/my-smoke
.venv/bin/dnastore benchmark configs/benchmark_badread_smoke.yaml \
  --out results/phase4/my-badread-smoke
.venv/bin/dnastore benchmark configs/benchmark_validation.yaml \
  --out results/phase4/my-validation
.venv/bin/dnastore benchmark configs/paper.yaml \
  --out results/phase4/my-study
```

The runner matches **exact total synthesized nucleotides**, including headers and
parity, across all arms. It records noiseless recovery, independent channel trials,
both decoder strengths, Wilson confidence intervals, empirical depth-search results,
continuous error metrics, PNG/SVG plots and provenance. Study mode requires at least
20 channel seeds per probe. The output directory must be new or empty.

See [`docs/BENCHMARK.md`](docs/BENCHMARK.md) for configuration, fairness, statistics and
runtime limits. Badread is the main read-level channel; IDS is the controlled reference.
Squigulator → Dorado remains a miscalibrated stress channel. Calibration caveats are in
[`docs/notes/channel_calibration.md`](docs/notes/channel_calibration.md).

## Project map and validation

- `dnastore/`: codecs, container, ECC, encoder/decoder, channels, screening and benchmark.
- `configs/`: encoder settings, primers, six-arm smoke and larger study presets.
- `scripts/`: calibration, characterization, raw-signal and clustering-scale experiments.
- `tests/`: unit, property, CLI, file-format and end-to-end recovery checks.
- `docs/notes/`: methodology, calibration, assumptions and historical handoffs.
- `results/`: local experiment output; large raw signal assets are optional.

```sh
.venv/bin/pytest -q
.venv/bin/python -m pip check
```

Current implementation, validation results and remaining research work are tracked in
[`docs/PROJECT_STATUS.md`](docs/PROJECT_STATUS.md). Historical notes preserve earlier
format versions and hypotheses; the latest DESIGN.md revisions and status file take
precedence.
