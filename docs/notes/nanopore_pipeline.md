# Raw-signal stress channel: squigulator → Dorado

Phase 3, 2026-10-02. Code: `dnastore/nanopore.py`. Setup: `scripts/setup_nanopore_tools.sh`. Benchmark: `results/phase3/nanopore/bench_nanopore.py`. Raw results: `results/phase3/nanopore/` (`throughput.json`, `accuracy.csv`, `noise_decomposition.*`, `error_context.json`).

## Versions (`tools/VERSIONS.json`)
| Tool | Version | Archive sha256 |
|---|---|---|
| squigulator | 0.5.0 (prebuilt x86_64) | `7197920bc50b2c10638e61d47f85c9a1c90dc62d16125a096d51e995db803a7d` |
| Dorado | 2.1.2+8b8fc5d (latest stable on GitHub, 2026-08-26) | `f4ed83acfb75cf07ffe8a0fc78e26828fc911fcfc8177920be6104e1d0e02485` |
| Models | `dna_r10.4.1_e8.2_400bps_hac@v6.0.0`, `fast@v5.2.0`, `sup@v5.2.0` | via `dorado download` |
| blue-crab (BLOW5→POD5) | 0.6.0 (pod5 0.3.48, pyslow5 1.5.0), isolated in `tools/venv` | — |

- **GPU:** RTX 3050 Laptop, 4096 MiB, driver 580.178.04 (CUDA 13.0).
- **CDN host:** the ONT CDN is `cdn.oxfordnanoportal.com`. The host `cdn.oxfordnanopore.com` does not resolve (NXDOMAIN).

## Commands
```
squigulator -x dna-r10-prom reference.fa -o reads.blow5 --full-contigs --ont-friendly=yes --seed S -t 1 --prefix=no
blue-crab s2p reads.blow5 -o reads.pod5
dorado basecaller tools/models/dna_r10.4.1_e8.2_400bps_hac@v6.0.0 reads.pod5 --emit-fastq --no-trim --device cuda:0 > basecalls_hac.fastq
```
Every command and its stderr is appended to `<out>/commands.log`. Versions, seed, and timings go to `<out>/pipeline.json`.

## Chemistry and sampling match
- squigulator v0.5.0 `dna-r10-prom` uses **5 kHz** sampling (since v0.3.0; `--sample-rate 4000` gives the old 4 kHz). Its BLOW5 header reports SQK-LSK114 and FLO-PRO114M. Dwell is 13 samples per base, about 385 bases/s.
- The Dorado models are R10.4.1 e8.2 400 bps, which in the current lineage are 5 kHz models, so the two match.
- **squigulator warns that its "parameters and models for dna-r10-prom 5khz are still crude".** Its R10 pore model is derived from nanoporetech/kmer_models.
- Dorado does not read BLOW5, so the files are converted with blue-crab `s2p`. POD5 requires UUID read IDs, so squigulator runs with `--ont-friendly=yes`. That gives sequential UUIDs `00000000-0000-0000-0000-{i+1:012d}` in FASTA record order, and `truth.tsv` maps them back to (oligo, copy, strand).

## Strand and copy number
`--full-contigs` simulates every FASTA record once, on the + strand. Copy number is applied by writing one record per molecule copy. Orientation is a seeded per-copy choice (P(rc) = 0.5): the record holds the oligo or its reverse complement.

## Short-read handling
- squigulator's 9-mer model emits signal for len − (k−1) bases (a 44-nt contig gives 36 simulated bases), so the ends of a bare 200-nt oligo lose primer bases.
- Each record is therefore wrapped in constant flanks: `FLANK_5` = poly-T + the ONT Y-adapter top-strand motif `CCTGTACTTCGTTCAGTTACGTATTGCT` as listed in Porechop/Dorado trimming tables, and `FLANK_3` = its reverse complement. They are *adapter-like*, not verified to equal the exact SQK-LSK114 adapter, and carry no data.
- Dorado runs with `--no-trim`, so it does not cut our primers. Primer finding belongs to our decoder.

200 oligos, hac:

| Variant | Edit rate | Full length (edit ≤ 20%) |
|---|---|---|
| No flanks | 13.96% | 96% |
| Flanks (**default**) | 10.26% | 100% |
| squigulator `--prefix=yes` | 10.57% | 100% |

No reads were dropped in any variant.

## Throughput (flanks; distinct random F + 160 nt + rc(R) oligos, one molecule each)
| Model | Reads | squigulator | blue-crab | Dorado | Dorado reads/s | Samples/s | Peak VRAM |
|---|---|---|---|---|---|---|---|
| hac v6.0.0 | 2,000 | 2.1 s | 2.6 s | 9.4 s | 213 | 0.71 M | 2.4 GB |
| hac v6.0.0 | 10,000 | 6.5 s | 13.1 s | 12.7 s | 788 | 2.6 M | 2.4 GB |
| hac v6.0.0 | 50,000 | 30.0 s | 45.6 s | 35.2 s | 1,422 | 4.7 M | 3.6 GB |
| sup v5.2.0 | 2,000 | 1.9 s | 4.2 s | 14.9 s | 135 | 0.45 M | 2.5 GB |
| sup v5.2.0 | 10,000 | 7.5 s | 10.2 s | 55.3 s | 181 | 0.60 M | 2.7 GB |
| fast v5.2.0 | 2,000 | 1.6 s | 1.7 s | 10.5 s | 190 | 0.63 M | 0.9 GB |
| fast v5.2.0 | 10,000 | 5.7 s | 9.6 s | 11.9 s | 842 | 2.8 M | 0.9 GB |

- Dorado's fixed start-up of about 7–8 s dominates small runs. The end-to-end cost is about **111 s per 50k reads** with hac.
- **sup is feasible** on 4 GB, at roughly 8× slower Dorado time than hac at 10k reads.

## Accuracy (per read vs source oligo, edlib infix alignment of the oligo inside the read)
| Model | n | Mean edit rate | Sub | Ins | Del | Median identity | Full length |
|---|---|---|---|---|---|---|---|
| hac | 50,000 | 9.67% | 4.19% | 3.12% | 2.35% | 0.907 | 99.72% |
| sup | 10,000 | 8.72% | 3.85% | 2.84% | 2.03% | 0.917 | 99.69% |
| fast | 10,000 | 15.85% | 7.23% | 5.83% | 2.80% | 0.848 | 86.3% |

**Error decomposition (hac, 2000 oligos, squigulator noise switches):**

| squigulator setting | Edit rate | Full length |
|---|---|---|
| Ideal (no amplitude or time noise) | 30.5% | 4.7% |
| Ideal amplitude | 31.3% | — |
| Ideal time | 8.5% | — |
| Default | 9.7% | — |

- Noise-free amplitude is badly out of distribution for Dorado.
- Most of the ~9% floor is *not* squigulator's added noise. It is the mismatch between squigulator's k-mer-model signal and the real signal Dorado was trained on.
- These initial measurements describe the simulated channel. Subsequent comparison with real R10.4.1 data is in `channel_calibration.md` and the calibration update below.

**Does channel error track our objective?** Spearman ρ, per-oligo mean edit rate (hac, 10k oligos) vs oligo features (bootstrap 95% CI):

| Feature | ρ [95% CI] |
|---|---|
| R10 Hamming-neighbour confusability, mean (δ = 0.05) | 0.004 [−0.014, 0.023] |
| R10 Hamming-neighbour confusability, max | 0.002 [−0.017, 0.023] |
| Whritenour-style adjacent-level contrast (R10, fraction of steps ≤ τ = 0.5 std units) | **0.13 [0.11, 0.15]** |
| GC fraction | −0.015 [−0.036, 0.004] |

In this channel our Hamming-neighbour confusability score does **not** predict per-oligo error, while adjacent-level contrast does, weakly. This bears directly on the confusability objective in the beam, and on the user's request to optimise against R10.

## Limitations
1. **Partial circularity:** squigulator's R10 signal derives from ONT's k-mer level model (nanoporetech/kmer_models), the same family as our R10 confusability table. Dorado's learned error process is the independent part. Because the error floor comes from model mismatch, errors that depend on sequence context may still partly reflect the level table. The absence of correlation with our Hamming score shows that the objective is not trivially rewarded here.
2. squigulator's R10 5 kHz parameters are self-described as crude. Absolute error rates are channel-specific.
3. Only sequencing is simulated: no library prep, adapter ligation efficiency, chimeras, or PCR. Copy number is an input, so it can come from `dnastore.coverage`.
4. The flanks are an adapter-like stand-in. Real reads carry adapter, barcode, and partial-read effects that are not modelled.
5. Dorado 2.1.2 notes a non-determinism fix for the v6.0.0 hac model. Repeat runs should still be checked for bit-identical FASTQ before claiming exact reproducibility (not yet tested).

## Calibration update (2026-10-03)

The real-data calibration is documented in [`channel_calibration.md`](channel_calibration.md).
Against matched real R10.4.1 encoded intervals, squigulator → Dorado hac has 9.12% raw
per-base error versus 1.51% in real hac, and its systematic-site fraction at short run
lengths is about 300× higher. The consensus floor is therefore a simulator artefact, not
evidence that our codec arms fail on real nanopore data. Use calibrated Badread as the
read-level comparison channel and uniform IDS as the controlled channel; retain squigulator
as a stress test. Badread matches the raw error and run-length rates but has no systematic
errors between reads, unlike the small real component observed for homopolymers of length 5
or more. Calibration uses plasmid fragments rather than a 200-nt oligo pool, so this channel
choice remains an approximation.

## End-to-end decode through this channel (2026-10-02): systematic error floor
Script: `scripts/nanopore_e2e.py`. Diagnostics: `results/phase3/e2e_diag/` (smoke run, 800 B per arm, about 40 lognormal copies per oligo (σ = 0.27), hac v6.0.0).

**Decoder changes the run forced** (tests pass, 202/202):
1. **Primer search window.** `locate_and_orient` searched only len(primer) + 12 nt at each end, which rejected 100% of untrimmed (flanked) reads. The default margin is now 80 nt. Regression test: flanked reads on both strands, and no hits on 200 random 264-nt reads.
2. **Grouping.** Exact-index grouping fails at this error rate. Only **26%** of reads decode to their true index, 42% to an unused index, 1.5% to another *valid* index, 22% are invalid, and 7% are not oriented (consistent with 0.9^15 ≈ 0.2).
   - The new default is `grouping="cluster"`: `consensus.kmer_cluster` (k = 10 candidate lookup plus a bounded edlib check, read-to-read bound 0.30·L). The index is read from each cluster's candidate bodies, and sub-groups by exact read index are also tried, so a merged cluster can still yield several oligos.
   - The old path is still available as `grouping="index"`.
   - The bound is measured: same-oligo Dorado read pairs have median 0.156, p99 0.287, max 0.33. Distinct oligo bodies have min 0.36 (rotating ternary, median 0.425) and median 0.49 (whitened quaternary).
   - Clustering on the smoke run gives purity 1.00, and every oligo's reads fall into one cluster.

**Finding: consensus does not converge to the true sequence.** The residual errors are shared across reads, i.e. sequence-determined.

| Setting (30 random oligos × 30 copies) | Raw body error | Consensus error per nt | Exact consensus bodies |
|---|---|---|---|
| dna-r10-prom + hac v6.0.0 | 9.8% | 1.71% | 20% |
| dna-r10-prom + sup v5.2.0 | 9.2% | 1.50% | 20% |
| dna-r10-min + hac v6.0.0 | 9.9% | 1.54% | 23% |
| dna-r10-min + sup v5.2.0 | 9.2% | 1.54% | 7% |

- The same consensus code reaches ≥ 0.90 exact at 10% i.i.d. IDS and 20× (channel.md). A 30-read majority sharing a wrong base is not chance.
- 4 kHz Dorado models (≤ v4.1), which would match squigulator's older and better-validated 4 kHz R10 profile, are no longer offered (`dorado download --list-yaml` lists v4.2.0+ only, all 5 kHz).

Consensus floor by arm (all reads, about 40 per oligo, true grouping):

| Arm | Mean consensus edits per 160 nt | Exact |
|---|---|---|
| whiten 2-bit | 3.9 | 11% |
| fountain (screened, run ≤ 3) | 2.5 | 22% |
| goldman (run = 1) | 1.9 | 28% |

- Whitened arm, consensus error per nt by the true run length at the site: 1.9% (run 1), 2.1% (2), 2.8% (3), **6.5% (4), 25% (5)**.
- Error types: deletions > substitutions > insertions. Errors are spread evenly along the body.

**Consequence.** With per-oligo CRC-discard plus outer-RS erasure decoding (15% parity), **no arm recovers a file at any depth** in this channel: 72–89% of oligos carry at least one consensus error. File recovery and minimum depth for 100% recovery (primary metric b) are undefined here until an inner code, or error (not erasure) decoding, is added. Post-consensus per-oligo error (primary metric a) *is* measurable, and it differs by arm, in the direction of the constraints.

**Caveat (circularity, 0.3).** The floor is most plausibly squigulator's "crude" 5 kHz R10 model being mismatched to Dorado's training data. Its *context dependence* (homopolymers ≥ 4) may be a simulator artefact rather than a property of real R10.4.1 reads. Real R10.4.1 consensus accuracy is far higher, and has not been compared here.
