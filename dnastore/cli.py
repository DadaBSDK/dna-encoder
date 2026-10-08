"""Command-line interface (typer).

Commands: ``encode``, ``decode``, ``simulate``, ``primers-design``, ``primers-validate``.
``benchmark`` runs matched-budget, multi-seed experiments.
"""

from __future__ import annotations

import json
from pathlib import Path
from functools import wraps
from typing import Optional

import numpy as np
import typer
import yaml

from .config import load_config
from .decoder import decode_reads, read_fasta
from .encoder import encode_bytes, write_fasta
from .metrics import pool_metrics
from .primers import design_primer_pair, load_primers, save_primers, validate_primer_pair
from .screening import ScreeningConfig, screen_pool

app = typer.Typer(add_completion=False, help="Structure- and readout-aware DNA storage codec.")


def _user_errors(fn):
    @wraps(fn)
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (ValueError, OSError, yaml.YAMLError) as exc:
            typer.echo(f"Error: {exc}", err=True)
            raise typer.Exit(2) from exc
    return wrapped


def _summ(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=float)
    return {"min": float(a.min()), "mean": float(a.mean()), "max": float(a.max())} if a.size else {}


@app.command()
@_user_errors
def encode(
    input: Path = typer.Argument(..., exists=True, dir_okay=False, help="File to encode."),
    out_dir: Path = typer.Option(Path("results/encoded"), help="Output directory."),
    config: Optional[Path] = typer.Option(None, help="YAML config (defaults if omitted)."),
    codec: Optional[str] = typer.Option(None, help="Override codec name."),
    mfe: Optional[bool] = typer.Option(None, help="Compute ViennaRNA MFE during screening."),
    allow_violations: bool = typer.Option(False, help="Write the pool even if oligos violate hard constraints."),
) -> None:
    """Encode a file into an oligo pool (FASTA) plus a JSON manifest."""
    over: dict = {}
    if codec:
        over["codec"] = {"name": codec, "params": {}}
    if mfe is not None:
        over["compute_mfe"] = mfe
    cfg = load_config(config, over)
    primers = load_primers(cfg["primers"])
    res = encode_bytes(input.read_bytes(), primers, cfg)
    scfg = ScreeningConfig.from_dict(cfg["screening"])
    screens = screen_pool([s for _, s in res.oligos], primers, scfg, cfg["compute_mfe"], cfg.get("workers"))
    violations = [{"index": idx, "violations": s["violations"]} for (idx, _), s in zip(res.oligos, screens) if s["violations"]]
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{input.name}.{res.codec.name}"
    manifest = {
        "input": str(input),
        "config": cfg,
        "container_reason": res.info.reason,
        "metrics": pool_metrics(res),
        "screening": {
            "gc": _summ([s["gc"] for s in screens]),
            "max_homopolymer": _summ([s["max_hp"] for s in screens]),
            "mfe_min_both": _summ([min(s["mfe_fwd"], s["mfe_rc"]) for s in screens if s["mfe_fwd"] is not None]),
            "n_violating_oligos": len(violations),
            "frac_violating_oligos": len(violations) / max(1, len(screens)),
        },
        "violations": violations,
    }
    (out_dir / f"{stem}.manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    if violations and not allow_violations:
        typer.echo(f"{len(violations)} oligos violate hard constraints; pool NOT written (use --allow-violations).", err=True)
        raise typer.Exit(2)
    write_fasta(res.oligos, str(out_dir / f"{stem}.fasta"))
    m = manifest["metrics"]
    typer.echo(
        f"{m['n_oligos']} oligos ({m['n_header_oligos']} header, {m['n_data_oligos']} data+parity, {m['n_parity_rows']} parity rows); "
        f"code density {m['code_density']:.3f} b/nt; effective {m['effective_density']:.3f} b/nt; "
        f"{len(violations)} oligos with hard-constraint violations -> {out_dir / stem}.fasta"
    )


@app.command()
@_user_errors
def decode(
    reads: Path = typer.Argument(..., exists=True, dir_okay=False, help="FASTA/FASTQ reads."),
    out: Path = typer.Option(Path("results/decoded/recovered"), help="Output path; the extension comes from the header."),
    config: Optional[Path] = typer.Option(None, help="YAML config (primers, oligo_len)."),
    strength: str = typer.Option("repair", help="Decoder strength: erasure or repair."),
    grouping: str = typer.Option("cluster", help="Read grouping: cluster or index."),
) -> None:
    """Decode reads back into the original file."""
    cfg = load_config(config)
    res = decode_reads(read_fasta(str(reads)), load_primers(cfg["primers"]), cfg["oligo_len"],
                       workers=cfg.get("workers"), strength=strength, grouping=grouping)
    typer.echo(json.dumps(res.report, indent=2, default=str))
    if not res.ok:
        typer.echo("DECODE FAILED", err=True)
        raise typer.Exit(1)
    path = out.with_suffix(res.extension)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(res.data)
    typer.echo(f"recovered {len(res.data)} bytes -> {path} (CRC32 verified)")


@app.command("stream-encode")
@_user_errors
def stream_encode(
    input: Path = typer.Argument(..., exists=True, dir_okay=False, help="Any file: image, video, GIF or other bytes."),
    out: Path = typer.Option(..., help="New .dna archive path (never overwritten)."),
) -> None:
    """Encode any-size files as streaming DNA text, without error correction."""
    from .stream import encode_file
    typer.echo(json.dumps(encode_file(input, out), indent=2))


@app.command("stream-decode")
@_user_errors
def stream_decode(
    input: Path = typer.Argument(..., exists=True, dir_okay=False, help="DNASTORE .dna archive."),
    out: Path = typer.Option(..., help="Exact recovered filename including extension; must not exist."),
) -> None:
    """Recover exact bytes from a DNA archive; verify SHA-256 before publishing output."""
    from .stream import decode_file
    typer.echo(json.dumps(decode_file(input, out), indent=2))


@app.command()
@_user_errors
def simulate(
    pool: Path = typer.Argument(..., exists=True, dir_okay=False, help="FASTA oligo pool from `encode`."),
    out: Path = typer.Option(Path("results/simulated/reads.fasta"), help="Output reads (FASTA)."),
    channel: Optional[Path] = typer.Option(None, help="YAML with ChannelParams fields."),
    seed: int = typer.Option(0, help="Channel RNG seed."),
) -> None:
    """Simulate sequencing reads from an oligo pool (uniform or confusability-weighted IDS channel).

    The confusability-weighted mode is a CIRCULAR sanity check, not a nanopore model (DESIGN.md 0.3).
    The main read-level evaluation uses calibrated Badread; squigulator is a stress channel."""
    import yaml

    from .channel import ChannelParams, simulate as run_channel

    params = ChannelParams.from_dict(yaml.safe_load(channel.read_text()) if channel else {})
    oligos = read_fasta(str(pool))
    reads, truth = run_channel(oligos, params, np.random.default_rng(seed))
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for i, r in enumerate(reads):
            fh.write(f">read_{i} oligo={r.oligo} rev={int(r.reverse)}\n{r.seq}\n")
    typer.echo(f"{len(reads)} reads from {len(oligos)} oligos ({int(truth.dropped.sum())} dropped) -> {out}")


@app.command()
@_user_errors
def benchmark(
    config: Path = typer.Argument(..., exists=True, dir_okay=False, help="Benchmark YAML."),
    out: Path = typer.Option(Path("results/benchmark"), help="New or empty output directory."),
) -> None:
    """Run equal-total-nucleotide experiments with independent channel seeds and Wilson CIs."""
    from .benchmark import run_benchmark

    run_benchmark(config, out, progress=typer.echo)
    typer.echo(f"Benchmark complete -> {out}")


@app.command("primers-design")
def primers_design(seed: int = 1, out: Path = Path("configs/primers.yaml")) -> None:
    """Design a validated primer pair deterministically from ``seed``."""
    pair, report = design_primer_pair(seed)
    save_primers(pair, out, report)
    typer.echo(f"F={pair.forward} R={pair.reverse} -> {out}")


@app.command("primers-validate")
def primers_validate(path: Path = Path("configs/primers.yaml")) -> None:
    """Re-validate a primer file."""
    rep = validate_primer_pair(load_primers(path))
    typer.echo(json.dumps(rep, indent=2, default=str))
    if not rep["ok"]:
        raise typer.Exit(1)


@app.command("web")
@_user_errors
def web(
    port: int = typer.Option(8765, min=1, max=65535, help="Local browser port."),
    out_dir: Path = typer.Option(Path("results/workbench"), help="Saved uploads and experiment artifacts."),
) -> None:
    """Start the local DNA testing frontend at http://127.0.0.1:8765."""
    from .web import serve
    serve(port=port, root=out_dir)


if __name__ == "__main__":
    app()
