"""Sequence screening: GC content, homopolymers, secondary structure (MFE), primer cross-hits.

All functions operate on the *full* oligo (primers included) unless stated otherwise.

MFE
---
ViennaRNA with DNA parameters (``RNA.params_load_DNA_Mathews2004()``, verified on
ViennaRNA 2.7.2). The parameter load is global per process, so :func:`_ensure_dna_params`
runs it lazily once per process, which makes multiprocessing workers safe. DNA MFE is not
strand-symmetric (G·T wobble pairs vs. C·A, and dangles differ), so both strands are
folded (:func:`mfe_both`).

Primer cross-hits
-----------------
Reads come from both strands, so each of F, R, rc(F), rc(R) is slid over every window of
the forward-strand oligo string, using Hamming distance. The true primer sites are
excluded: ``F`` at position 0 and ``rc(R)`` at position ``L - len(R)``.

3'-anchored hits model mispriming, which is driven by complementarity at the primer's 3'
end. A primer P can prime wherever the forward-strand string contains P's 3' end
(``P[-a:]``, annealing to the reverse strand) or its reverse complement (annealing to the
forward strand). ``rc(P[-a:]) == rc(P)[:a]``, so for the reverse-complement queries the
3' end of the primer corresponds to the **first** ``a`` bases of the rc string. Queries:
``F[-a:]``, ``R[-a:]``, ``rc(F)[:a]``, ``rc(R)[:a]``. True sites excluded:
``F[-a:]`` at ``len(F) - a`` and ``rc(R)[:a]`` at ``L - len(R)``.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, fields
from typing import Iterable

import numpy as np
import RNA

from dnastore.dnautil import revcomp
from dnastore.primers import PrimerPair


@dataclass(frozen=True)
class ScreeningConfig:
    """Constraint thresholds. ``hard`` lists the constraints that count as violations."""

    gc_global: tuple[float, float] = (0.40, 0.60)
    gc_window_size: int = 20
    gc_window: tuple[float, float] = (0.25, 0.75)
    max_homopolymer: int = 3
    mfe_temperature: float = 37.0
    mfe_threshold: float | None = None
    primer_max_mismatch: int = 3
    primer_anchor_len: int = 8
    primer_anchor_mismatch: int = 0
    hard: tuple[str, ...] = ("homopolymer", "gc_global", "gc_window", "primer")

    @classmethod
    def from_dict(cls, d: dict | None) -> "ScreeningConfig":
        """Build from a plain dict (YAML). Unknown keys raise; lists become tuples."""
        d = dict(d or {})
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ValueError(f"unknown screening keys: {sorted(unknown)}")
        for k in ("gc_global", "gc_window", "hard"):
            if k in d and d[k] is not None:
                d[k] = tuple(d[k])
        return cls(**d)


# --------------------------------------------------------------------------- simple metrics


def gc_fraction(seq: str) -> float:
    """Fraction of G/C bases (0.0 for an empty sequence)."""
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def gc_window_extremes(seq: str, w: int) -> tuple[float, float]:
    """(min, max) GC fraction over all length-``w`` windows; whole sequence if shorter."""
    if not seq:
        return (0.0, 0.0)
    if len(seq) <= w:
        g = gc_fraction(seq)
        return (g, g)
    arr = np.frombuffer(seq.encode(), np.uint8)
    isgc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.int32)
    cs = np.concatenate([[0], np.cumsum(isgc)])
    win = (cs[w:] - cs[:-w]) / w
    return (float(win.min()), float(win.max()))


def max_homopolymer(seq: str) -> int:
    """Length of the longest run of identical bases."""
    if not seq:
        return 0
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        if run > best:
            best = run
    return best


# --------------------------------------------------------------------------- MFE

_DNA_PARAMS_LOADED = False


def _ensure_dna_params() -> None:
    global _DNA_PARAMS_LOADED
    if not _DNA_PARAMS_LOADED:
        RNA.params_load_DNA_Mathews2004()
        _DNA_PARAMS_LOADED = True


def mfe(seq: str, temperature: float = 37.0) -> float:
    """Minimum free energy (kcal/mol) of ``seq`` as single-stranded DNA."""
    _ensure_dna_params()
    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(seq, md)
    _, e = fc.mfe()
    return float(e)


def mfe_both(seq: str, temperature: float = 37.0) -> tuple[float, float]:
    """MFE of the sequence and of its reverse complement."""
    return mfe(seq, temperature), mfe(revcomp(seq), temperature)


# --------------------------------------------------------------------------- primer hits

_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i


def _enc(seq: str) -> np.ndarray:
    return _CODE[np.frombuffer(seq.encode(), np.uint8)]


def _hamming_scan(text: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Hamming distance of ``query`` against every window of ``text``."""
    if len(query) > len(text):
        return np.zeros(0, dtype=np.int32)
    win = np.lib.stride_tricks.sliding_window_view(text, len(query))
    return (win != query[None, :]).sum(axis=1)


def primer_hits(oligo: str, primers: PrimerPair, max_mm: int, anchor_len: int, anchor_mm: int) -> dict:
    """Spurious full-primer and 3'-anchored primer matches in both orientations.

    Returns ``{"full": n, "anchor": n, "full_positions": [(query, pos, mm)], "anchor_positions": [...]}``.
    """
    text = _enc(oligo)
    L = len(oligo)
    o = primers.orientations()
    lf, lr = len(primers.forward), len(primers.reverse)
    full_true = {("F", 0), ("rcR", L - lr)}
    full_pos = []
    for name, q in o.items():
        d = _hamming_scan(text, _enc(q))
        for pos in np.nonzero(d <= max_mm)[0]:
            if (name, int(pos)) not in full_true:
                full_pos.append((name, int(pos), int(d[pos])))
    a = anchor_len
    anchors = {
        "F3": primers.forward[-a:],
        "R3": primers.reverse[-a:],
        "rcF3": o["rcF"][:a],
        "rcR3": o["rcR"][:a],
    }
    anchor_true = {("F3", lf - a), ("rcR3", L - lr)}
    anchor_pos = []
    for name, q in anchors.items():
        d = _hamming_scan(text, _enc(q))
        for pos in np.nonzero(d <= anchor_mm)[0]:
            if (name, int(pos)) not in anchor_true:
                anchor_pos.append((name, int(pos), int(d[pos])))
    return {
        "full": len(full_pos),
        "anchor": len(anchor_pos),
        "full_positions": full_pos,
        "anchor_positions": anchor_pos,
    }


# --------------------------------------------------------------------------- confusability
# Confusability and adjacent-level contrast are implemented in dnastore.kmer;
# scripts/confusability_report.py provides the model and delta-sweep reports.


# --------------------------------------------------------------------------- screening


def screen_oligo(oligo: str, primers: PrimerPair, cfg: ScreeningConfig, compute_mfe: bool) -> dict:
    """All screening metrics for one full oligo, plus the list of hard-constraint violations."""
    gmin, gmax = gc_window_extremes(oligo, cfg.gc_window_size)
    gc = gc_fraction(oligo)
    hp = max_homopolymer(oligo)
    ph = primer_hits(oligo, primers, cfg.primer_max_mismatch, cfg.primer_anchor_len, cfg.primer_anchor_mismatch)
    mf = mr = None
    if compute_mfe:
        mf, mr = mfe_both(oligo, cfg.mfe_temperature)
    checks = {
        "homopolymer": hp <= cfg.max_homopolymer,
        "gc_global": cfg.gc_global[0] <= gc <= cfg.gc_global[1],
        "gc_window": cfg.gc_window[0] <= gmin and gmax <= cfg.gc_window[1],
        "primer": ph["full"] == 0,
        "primer_anchor": ph["anchor"] == 0,
    }
    if compute_mfe and cfg.mfe_threshold is not None:
        checks["mfe"] = min(mf, mr) >= cfg.mfe_threshold
    violations = [name for name in cfg.hard if name in checks and not checks[name]]
    return {
        "gc": gc,
        "gc_win_min": gmin,
        "gc_win_max": gmax,
        "max_hp": hp,
        "primer_hits": ph["full"],
        "primer_anchor_hits": ph["anchor"],
        "mfe_fwd": mf,
        "mfe_rc": mr,
        "violations": violations,
    }


def _screen_star(args: tuple) -> dict:
    return screen_oligo(*args)


def screen_pool(
    oligos: Iterable[str],
    primers: PrimerPair,
    cfg: ScreeningConfig,
    compute_mfe: bool,
    workers: int | None = None,
) -> list[dict]:
    """Screen many oligos. Uses a process pool when MFE is computed (the bottleneck)."""
    oligos = list(oligos)
    jobs = [(o, primers, cfg, compute_mfe) for o in oligos]
    if not compute_mfe or len(oligos) < 8 or workers == 1:
        return [_screen_star(j) for j in jobs]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        chunk = max(1, len(jobs) // ((workers or 8) * 4))
        return list(ex.map(_screen_star, jobs, chunksize=chunk))
