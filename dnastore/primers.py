"""Primer pair design and validation.

Convention: the physical oligo is ``forward + body + revcomp(reverse)``. ``reverse`` is the
reverse primer as ordered (5'->3'). It anneals to the forward strand at the 3' tail,
whose sequence is :attr:`PrimerPair.tail`.

Thermodynamics come from primer3-py (SantaLucia nearest-neighbour). ``calc_tm`` is called
with the primer3 defaults, which are recorded in :data:`PRIMER3_CONDITIONS`. Hairpin and
dimer ΔG are in cal/mol at 37 °C (primer3 default).

Cross-similarity check (kept simple, documented): for every pair of the four orientation
strings {F, R, rc(F), rc(R)}, including each string against its own reverse complement
(self-complementarity), slide one over the other at every shift with an overlap of at
least ``min_overlap`` nt. Record the mismatch fraction in the overlap. The pair passes if
(a) the minimum mismatch fraction over all shifts is ``>= min_mismatch_frac``, and (b) at
shift 0 (full-length alignment) the Hamming distance is ``>= min_hamming_full``. Pairs
that are reverse complements of each other by definition (X vs rc(X) at shift 0, i.e.
F/rcF and R/rcR) are compared at non-zero shifts only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import primer3
import yaml

from dnastore.dnautil import BASES, revcomp

PRIMER3_CONDITIONS = {
    "mv_conc_mM": 50.0,
    "dv_conc_mM": 1.5,
    "dntp_conc_mM": 0.6,
    "dna_conc_nM": 50.0,
    "tm_method": "santalucia",
    "salt_corrections_method": "santalucia",
    "dG_temperature_C": 37.0,
}


@dataclass(frozen=True)
class PrimerPair:
    """Forward and reverse primers (both 5'->3' as ordered)."""

    forward: str
    reverse: str

    @property
    def tail(self) -> str:
        """Sequence appended to the 3' end of the forward strand: ``revcomp(reverse)``."""
        return revcomp(self.reverse)

    def orientations(self) -> dict[str, str]:
        """The four strings that can appear in reads: F, R, rc(F), rc(R)."""
        return {
            "F": self.forward,
            "R": self.reverse,
            "rcF": revcomp(self.forward),
            "rcR": revcomp(self.reverse),
        }


@dataclass(frozen=True)
class PrimerConstraints:
    """Acceptance criteria for a primer pair."""

    length: int = 20
    gc_min: float = 0.45
    gc_max: float = 0.55
    tm_min: float = 58.0
    tm_max: float = 62.0
    max_tm_diff: float = 2.0
    max_homopolymer: int = 3
    clamp_window: int = 5
    clamp_gc_min: int = 1
    clamp_gc_max: int = 3
    min_hairpin_dg: float = -2000.0
    min_dimer_dg: float = -6000.0
    min_overlap: int = 12
    min_mismatch_frac: float = 0.4
    min_hamming_full: int = 10


def _max_run(seq: str) -> int:
    best = run = 1 if seq else 0
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def _shift_mismatch(a: str, b: str, min_overlap: int, skip_zero: bool) -> tuple[float, int | None]:
    """Minimum mismatch fraction over shifts with overlap >= ``min_overlap``; also return the
    Hamming distance at shift 0 (None if lengths differ or ``skip_zero``)."""
    x = np.frombuffer(a.encode(), np.uint8)
    y = np.frombuffer(b.encode(), np.uint8)
    best = 1.0
    h0 = None
    for s in range(-(len(y) - min_overlap), len(x) - min_overlap + 1):
        if s == 0 and skip_zero:
            continue
        xs = x[max(s, 0) : max(s, 0) + len(y) - max(-s, 0)]
        ys = y[max(-s, 0) : max(-s, 0) + len(xs)]
        n = min(len(xs), len(ys))
        if n < min_overlap:
            continue
        mm = int((xs[:n] != ys[:n]).sum())
        best = min(best, mm / n)
        if s == 0 and len(x) == len(y):
            h0 = mm
    return best, h0


def _tm(seq: str) -> float:
    return float(primer3.calc_tm(seq))


def _single_checks(seq: str, c: PrimerConstraints) -> dict:
    tail = seq[-c.clamp_window :]
    clamp = tail.count("G") + tail.count("C")
    return {
        "seq": seq,
        "gc": _gc(seq),
        "tm": _tm(seq),
        "max_homopolymer": _max_run(seq),
        "clamp_gc": clamp,
        "hairpin_dg": float(primer3.calc_hairpin(seq).dg),
        "homodimer_dg": float(primer3.calc_homodimer(seq).dg),
    }


def _single_failures(name: str, m: dict, c: PrimerConstraints) -> list[str]:
    f = []
    if len(m["seq"]) != c.length:
        f.append(f"{name}: length {len(m['seq'])}")
    if not c.gc_min <= m["gc"] <= c.gc_max:
        f.append(f"{name}: gc {m['gc']:.2f}")
    if not c.tm_min <= m["tm"] <= c.tm_max:
        f.append(f"{name}: tm {m['tm']:.1f}")
    if m["max_homopolymer"] > c.max_homopolymer:
        f.append(f"{name}: homopolymer {m['max_homopolymer']}")
    if not c.clamp_gc_min <= m["clamp_gc"] <= c.clamp_gc_max:
        f.append(f"{name}: 3' clamp GC {m['clamp_gc']}")
    if m["hairpin_dg"] <= c.min_hairpin_dg:
        f.append(f"{name}: hairpin dG {m['hairpin_dg']:.0f}")
    if m["homodimer_dg"] <= c.min_dimer_dg:
        f.append(f"{name}: homodimer dG {m['homodimer_dg']:.0f}")
    return f


def validate_primer_pair(pair: PrimerPair, c: PrimerConstraints | None = None) -> dict:
    """Compute all metrics for ``pair``. Returns a dict with ``ok`` and ``failures``."""
    c = c or PrimerConstraints()
    fm = _single_checks(pair.forward, c)
    rm = _single_checks(pair.reverse, c)
    failures = _single_failures("F", fm, c) + _single_failures("R", rm, c)
    tm_diff = abs(fm["tm"] - rm["tm"])
    if tm_diff > c.max_tm_diff:
        failures.append(f"|dTm| {tm_diff:.1f}")
    hetero = float(primer3.calc_heterodimer(pair.forward, pair.reverse).dg)
    if hetero <= c.min_dimer_dg:
        failures.append(f"heterodimer dG {hetero:.0f}")
    # Cross-similarity among orientation strings.
    o = pair.orientations()
    names = list(o)
    cross = {}
    for i in range(len(names)):
        for j in range(i, len(names)):
            a, b = names[i], names[j]
            if a == b:
                continue
            skip0 = {a, b} in ({"F", "rcF"}, {"R", "rcR"})
            frac, h0 = _shift_mismatch(o[a], o[b], c.min_overlap, skip0)
            cross[f"{a}~{b}"] = {"min_mismatch_frac": round(frac, 3), "hamming_shift0": h0}
            if frac < c.min_mismatch_frac:
                failures.append(f"cross {a}~{b}: mismatch frac {frac:.2f}")
            if h0 is not None and h0 < c.min_hamming_full:
                failures.append(f"cross {a}~{b}: full Hamming {h0}")
    return {
        "ok": not failures,
        "failures": failures,
        "forward": fm,
        "reverse": rm,
        "tm_diff": tm_diff,
        "heterodimer_dg": hetero,
        "cross": cross,
        "constraints": asdict(c),
        "primer3_conditions": PRIMER3_CONDITIONS,
    }


def _random_candidate(rng: np.random.Generator, c: PrimerConstraints) -> str | None:
    seq = "".join(BASES[i] for i in rng.integers(0, 4, c.length))
    m = {"seq": seq, "gc": _gc(seq), "max_homopolymer": _max_run(seq)}
    if not c.gc_min <= m["gc"] <= c.gc_max or m["max_homopolymer"] > c.max_homopolymer:
        return None
    tail = seq[-c.clamp_window :]
    if not c.clamp_gc_min <= tail.count("G") + tail.count("C") <= c.clamp_gc_max:
        return None
    full = _single_checks(seq, c)
    return seq if not _single_failures("X", full, c) else None


def design_primer_pair(seed: int, c: PrimerConstraints | None = None, max_tries: int = 200_000) -> tuple[PrimerPair, dict]:
    """Deterministic rejection sampling of a primer pair satisfying ``c``."""
    c = c or PrimerConstraints()
    rng = np.random.Generator(np.random.PCG64(seed))
    singles: list[str] = []
    for tries in range(1, max_tries + 1):
        s = _random_candidate(rng, c)
        if s is None:
            continue
        for f in singles:
            pair = PrimerPair(f, s)
            rep = validate_primer_pair(pair, c)
            if rep["ok"]:
                rep["design_seed"] = seed
                rep["candidates_sampled"] = tries
                return pair, rep
        singles.append(s)
    raise RuntimeError(f"no valid primer pair found in {max_tries} samples")


def save_primers(pair: PrimerPair, path: str | Path, report: dict) -> None:
    """Write primers and their validation report to YAML."""
    doc = {
        "forward": pair.forward,
        "reverse": pair.reverse,
        "design_seed": report.get("design_seed"),
        "validation": report,
    }
    Path(path).write_text(yaml.safe_dump(doc, sort_keys=False))


def load_primers(path: str | Path) -> PrimerPair:
    """Read a primer pair from YAML written by :func:`save_primers` (or any YAML with
    ``forward`` and ``reverse`` keys)."""
    doc = yaml.safe_load(Path(path).read_text())
    return PrimerPair(str(doc["forward"]).upper(), str(doc["reverse"]).upper())
