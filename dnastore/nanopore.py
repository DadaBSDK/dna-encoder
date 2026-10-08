"""Main evaluation channel: squigulator raw-signal simulation -> Dorado basecalling.

Pipeline (docs/notes/nanopore_pipeline.md)::

    oligos x copies --(seeded strand choice, optional flanks)--> reference FASTA
      --squigulator -x dna-r10-prom --full-contigs --ont-friendly=yes--> BLOW5 (5 kHz)
      --blue-crab s2p--> POD5
      --dorado basecaller <R10.4.1 e8.2 400bps 5kHz model> --emit-fastq --no-trim--> FASTQ

* **Strand.** squigulator ``--full-contigs`` simulates every FASTA record on the + strand,
  so orientation is chosen here: each copy is written as the oligo or its reverse
  complement with probability ``rc_frac`` (seeded). The truth table records the choice.
* **Read IDs.** ``--ont-friendly=yes`` gives sequential UUIDs
  ``00000000-0000-0000-0000-{i+1:012d}`` in FASTA record order (needed because POD5 requires
  UUIDs). The truth table maps them back to (oligo index, copy, strand).
* **Short reads.** squigulator's 9-mer model emits signal for ``len - (k-1)`` bases, and
  basecallers are least accurate at read ends. Each record can therefore be wrapped in
  fixed flanks (``FLANK_5`` / ``FLANK_3``, constant, not part of the oligo), and/or use
  squigulator's adapter prefix (``prefix=True``). The flanks stand in for the sequencing
  adapter / ligation context that a real library adds. They are not stored data, and the
  decoder finds primers inside the read anyway.
* **Determinism.** squigulator output is bit-identical for a fixed ``--seed`` only with
  ``-t 1`` (with 8 threads the per-read RNG streams depend on scheduling; verified), so
  ``threads=1`` is the default. Dorado determinism is checked separately (see notes).
* **Circularity.** squigulator signal is generated from ONT k-mer level models (R10 from
  nanoporetech/kmer_models), the same family as our confusability table. Only Dorado's error
  process is independent of our objective.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .dnautil import revcomp

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"
SQUIGULATOR = TOOLS / "squigulator-v0.5.0" / "squigulator"
BLUE_CRAB = TOOLS / "venv" / "bin" / "blue-crab"
MODELS_DIR = TOOLS / "models"
PROFILE = "dna-r10-prom"  # 5 kHz, LSK114, FLO-PRO114M (matches Dorado e8.2 400bps 5kHz models)

# Constant flanks (no relation to data). FLANK_5 mimics a ligation-adapter context;
# FLANK_3 covers the (k-1) bases the 9-mer signal model cannot emit at the 3' end.
# FLANK_5: poly-T leader + the ONT ligation Y-adapter top-strand motif as listed by
# Porechop/Dorado trimming tables ("CCTGTACTTCGTTCAGTTACGTATTGCT"). It is adapter-*like*:
# not verified to equal the exact SQK-LSK114 adapter, and used only as a constant context.
FLANK_5 = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"
FLANK_3 = "AGCAATACGTAACTGAACGAAGTACAGG"  # revcomp of the motif


def dorado_bin() -> Path | None:
    hits = sorted(glob.glob(str(TOOLS / "dorado-*" / "bin" / "dorado")))
    return Path(hits[-1]) if hits else None


def model_path(kind: str = "hac") -> Path | None:
    """Path of a downloaded R10.4.1 e8.2 400bps model of the given kind (fast/hac/sup)."""
    hits = sorted(p for p in glob.glob(str(MODELS_DIR / f"dna_r10.4.1_e8.2_400bps_{kind}@*")) if os.path.isdir(p))
    return Path(hits[-1]) if hits else None


def tools_available(kind: str = "hac") -> bool:
    return SQUIGULATOR.exists() and BLUE_CRAB.exists() and dorado_bin() is not None and model_path(kind) is not None


def _run(cmd: list[str], log: Path, stdout_path: Path | None = None) -> float:
    """Run a command, append it and its stderr to ``log``, return wall time (s)."""
    t0 = time.perf_counter()
    with open(log, "a") as lf:
        lf.write("$ " + " ".join(map(str, cmd)) + "\n")
        lf.flush()
        out = open(stdout_path, "w") if stdout_path else subprocess.DEVNULL
        try:
            proc = subprocess.run(list(map(str, cmd)), stdout=out, stderr=lf, check=False)
        finally:
            if stdout_path:
                out.close()
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(map(str, cmd))}; see {log}")
    return time.perf_counter() - t0


def tool_versions() -> dict[str, str]:
    """Versions of every external tool (recorded next to outputs)."""
    v: dict[str, str] = {}
    try:
        v["squigulator"] = subprocess.run([str(SQUIGULATOR), "--version"], capture_output=True, text=True).stdout.strip()
    except OSError:
        v["squigulator"] = "missing"
    d = dorado_bin()
    if d:
        r = subprocess.run([str(d), "--version"], capture_output=True, text=True)
        lines = [l for l in (r.stdout + r.stderr).splitlines() if l.strip() and not l.startswith("[")]
        v["dorado"] = lines[0].strip() if lines else "?"
    try:
        r = subprocess.run([str(BLUE_CRAB), "--version"], capture_output=True, text=True)
        v["blue_crab"] = (r.stdout + r.stderr).strip().splitlines()[0]
    except OSError:
        v["blue_crab"] = "missing"
    for kind in ("fast", "hac", "sup"):
        m = model_path(kind)
        if m:
            v[f"model_{kind}"] = m.name
    return v


def write_reference(oligos: Sequence[str], copies: Sequence[int], path: Path, seed: int, rc_frac: float = 0.5,
                    flanks: bool = True) -> list[dict[str, Any]]:
    """Write one FASTA record per molecule copy. Returns the truth table (record order)."""
    rng = np.random.default_rng(seed)
    truth = []
    with open(path, "w") as fh:
        for oi, (seq, n) in enumerate(zip(oligos, copies)):
            for c in range(int(n)):
                rev = bool(rng.random() < rc_frac)
                s = revcomp(seq) if rev else seq
                if flanks:
                    s = FLANK_5 + s + FLANK_3
                i = len(truth)
                fh.write(f">m{i}\n{s}\n")
                truth.append({"read_id": f"00000000-0000-0000-0000-{i + 1:012d}", "oligo": oi, "copy": c, "reverse": rev})
    return truth


def simulate_reads(oligos: Sequence[str], copies: Sequence[int], out_dir: str | Path, seed: int, profile: str = PROFILE,
                   threads: int = 1, flanks: bool = True, prefix: bool = False, rc_frac: float = 0.5,
                   extra_args: Sequence[str] = ()) -> dict[str, Any]:
    """Simulate raw signal for every copy of every oligo. Returns paths, truth, and timing."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log = out / "commands.log"
    ref = out / "reference.fa"
    truth = write_reference(oligos, copies, ref, seed, rc_frac, flanks)
    with open(out / "truth.tsv", "w") as fh:
        fh.write("read_id\toligo\tcopy\treverse\n")
        for t in truth:
            fh.write(f"{t['read_id']}\t{t['oligo']}\t{t['copy']}\t{int(t['reverse'])}\n")
    blow5 = out / "reads.blow5"
    if not truth:
        raise ValueError("no molecules to simulate")
    cmd = [SQUIGULATOR, "-x", profile, str(ref), "-o", blow5, "--full-contigs", "--ont-friendly=yes",
           "--seed", str(int(seed) or 1), "-t", str(threads), f"--prefix={'yes' if prefix else 'no'}", *extra_args]
    t = _run(cmd, log)
    return {"blow5": blow5, "truth": truth, "reference": ref, "t_squigulator_s": t}


def to_pod5(blow5: Path, out_dir: Path) -> tuple[Path, float]:
    pod5 = Path(out_dir) / "reads.pod5"
    if pod5.exists():
        pod5.unlink()
    t = _run([BLUE_CRAB, "s2p", str(blow5), "-o", str(pod5)], Path(out_dir) / "commands.log")
    return pod5, t


def basecall(signal_path: str | Path, out_dir: str | Path, model: str = "hac", device: str = "cuda:0",
             extra_args: Sequence[str] = ()) -> tuple[Path, float]:
    """Basecall a POD5 file (or directory) with Dorado; returns (fastq path, seconds).

    ``--no-trim`` keeps every basecalled base: our primers are not ONT kit primers, and
    trimming decisions belong to our decoder.
    """
    d = dorado_bin()
    m = model_path(model)
    if d is None or m is None:
        raise RuntimeError(f"dorado or the {model} model is not installed (run scripts/setup_nanopore_tools.sh)")
    out = Path(out_dir)
    fq = out / f"basecalls_{model}.fastq"
    cmd = [d, "basecaller", str(m), str(signal_path), "--emit-fastq", "--no-trim", "--device", device, *extra_args]
    t = _run(cmd, out / "commands.log", stdout_path=fq)
    return fq, t


def read_fastq(path: str | Path) -> list[tuple[str, str, str]]:
    """(read_id, sequence, quality) for every record; the read id is the first header token."""
    recs = []
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s, _, q = fh.readline().strip(), fh.readline(), fh.readline().strip()
            recs.append((h[1:].split()[0], s, q))
    return recs


def run_pipeline(oligos: Sequence[str], copies: Sequence[int], out_dir: str | Path, seed: int, model: str = "hac",
                 device: str = "cuda:0", flanks: bool = True, prefix: bool = False, keep_signal: bool = False,
                 threads: int = 1) -> tuple[Path, dict[str, Any]]:
    """Oligos + copy numbers -> basecalled FASTQ. Writes versions and timings next to the output."""
    out = Path(out_dir)
    sim = simulate_reads(oligos, copies, out, seed, flanks=flanks, prefix=prefix, threads=threads)
    pod5, t_conv = to_pod5(sim["blow5"], out)
    fq, t_bc = basecall(pod5, out, model=model, device=device)
    n_signal = sum(int(c) for c in copies)
    meta = {
        "versions": tool_versions(), "profile": PROFILE, "model": model, "seed": seed, "flanks": flanks,
        "prefix": prefix, "n_molecules": n_signal, "t_squigulator_s": sim["t_squigulator_s"],
        "t_s2p_s": t_conv, "t_dorado_s": t_bc,
    }
    (out / "pipeline.json").write_text(json.dumps(meta, indent=2, default=str))
    if not keep_signal:
        for p in (sim["blow5"], pod5):
            if Path(p).exists():
                Path(p).unlink()
    return fq, meta


def gpu_memory_used_mib() -> int | None:
    if shutil.which("nvidia-smi") is None:
        return None
    r = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"], capture_output=True, text=True)
    try:
        return int(r.stdout.strip().splitlines()[0])
    except (ValueError, IndexError):
        return None


# --------------------------------------------------------------------------- Badread channel
BADREAD = Path(sys.executable).parent / "badread"
# Calibrated to real R10.4.1 reads of synthetic encoded DNA (Chen et al. 2025, Zenodo
# 10.5281/zenodo.16883332) basecalled with Dorado hac v6.0.0: per-read identity mean 98.48%,
# sd 1.17 (docs/notes/channel_calibration.md).
BADREAD_IDENTITY = "98.48,99.9,1.2"


def badread_reads(oligos: Sequence[str], copies: Sequence[float], out_dir: str | Path, seed: int,
                  identity: str = BADREAD_IDENTITY, model: str = "nanopore2023",
                  flanks: bool = True) -> tuple[list[tuple[str, str]], list[dict[str, Any]]]:
    """Read-level empirical channel: Badread (v0.4.1) error + qscore models for R10.4.1.

    Each oligo (with the constant flanks) is a reference contig with relative depth
    ``copies[i]`` (Badread ``depth=`` header), and the expected total read count is
    ``sum(copies)``. Fragment lengths far exceed the contigs, so every read spans its whole
    molecule. No junk, random, chimeric, glitch, or adapter sequence is added. Strand is
    random per read (Badread). Returns ``(reads, truth)``, where truth rows are
    ``{read_id, oligo, reverse}``.

    Badread draws each read's errors independently from a context (k-mer) error model, so it
    has **no read-to-read systematic component**. The calibration found such a component in
    real reads only at homopolymer runs >= 5 (about 0.4-0.7% of positions).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ref = out / "badread_ref.fa"
    w = np.asarray(copies, dtype=float)
    keep = [i for i in range(len(oligos)) if w[i] > 0]
    total_bases = 0.0
    with open(ref, "w") as fh:
        for i in keep:
            s = (FLANK_5 + oligos[i] + FLANK_3) if flanks else oligos[i]
            fh.write(f">o{i} depth={w[i]:.6g}\n{s}\n")
            total_bases += w[i] * len(s)
    fq = out / "badread.fastq"
    cmd = [str(BADREAD), "simulate", "--reference", str(ref), "--quantity", str(int(round(total_bases))),
           "--error_model", model, "--qscore_model", model, "--identity", identity, "--length", "5000,500",
           "--junk_reads", "0", "--random_reads", "0", "--chimeras", "0", "--glitches", "0,0,0",
           "--start_adapter_seq", "", "--end_adapter_seq", "", "--seed", str(int(seed))]
    log = out / "commands.log"
    with open(fq, "w") as fo, open(log, "a") as lg:
        lg.write(" ".join(cmd) + "\n")
        subprocess.run(cmd, stdout=fo, stderr=lg, check=True)
    reads, truth = [], []
    for rid, seq, _ in read_fastq_full(fq):
        name, rest = rid
        contig, strand, _ = rest.split(",", 2)
        reads.append((name, seq))
        truth.append({"read_id": name, "oligo": int(contig[1:]), "reverse": strand.startswith("-")})
    return reads, truth


def read_fastq_full(path: str | Path) -> list[tuple[tuple[str, str], str, str]]:
    """Like :func:`read_fastq` but keeps the first two header tokens: ((id, second), seq, qual)."""
    recs = []
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s, _, q = fh.readline().strip(), fh.readline(), fh.readline().strip()
            tok = h[1:].split()
            recs.append(((tok[0], tok[1] if len(tok) > 1 else ""), s, q))
    return recs
