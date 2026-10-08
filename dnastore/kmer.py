"""Nanopore k-mer level models and the confusability / adjacent-contrast metrics.

Models (downloaded by ``scripts/fetch_kmer_models.py``; provenance and SHA-256 in
``data/kmer_models/PROVENANCE.json``; nanoporetech/kmer_models @ 4e56dae):

* ``r10.4.1_9mer``: ``dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt``. Expected level only,
  in **standardised units** (mean ~0, sd ~1). No per-k-mer std, no pA scale.
* ``r9.4.1_6mer``: ``legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model``.
  ``level_mean`` and ``level_stdv`` in **pA**. ONT names the directory "r9.4"; it is the
  model commonly used for R9.4.1 flow cells, and we label it r9.4.1 for that reason.

k-mer orientation assumption: k-mers are read 5'->3' along the strand that is basecalled.
Both strands of a dsDNA oligo are read, so every oligo-level metric is computed on the
sequence and on its reverse complement.

k-mer codes: A=0, C=1, G=2, T=3, most-significant base first, so
``code(x) = sum(b_i * 4**(k-1-i))``. Rolling update along a sequence::

    code = ((code << 2) | b) & (4**k - 1)

Reverse-complement strand lookup for incremental scoring: ``rc_code[c]`` is the code of
the reverse complement of k-mer ``c``. The k-mers of revcomp(seq) are exactly the reverse
complements of the k-mers of seq (in reverse order), so a sum or max over the RC strand
is ``table[rc_code[codes]]`` with no need to build the RC string.

Confusability (DESIGN.md §2.3, Part 0.6)::

    conf(x) = |{ y : d_H(x, y) = 1, |l(x) - l(y)| < delta(x, y) }| / (3k)   in [0, 1]

* R10.4.1: ``delta`` is a constant in standardised units (default 0.05).
* R9.4.1: ``delta(x, y) = c * sqrt((sd(x)**2 + sd(y)**2) / 2)`` in pA (default c = 1). This
  pooled (RMS) std makes the relation symmetric, so y confuses x iff x confuses y, and it
  is the usual scale for separating two Gaussians with unequal widths.

Adjacent contrast (Whritenour, Civelek & Farnoud 2025, Sci. Rep.,
doi:10.1038/s41598-025-08531-z): for consecutive overlapping k-mers (shift 1),
``|mu(x_i) - mu(x_{i+1})|``. Their constraint requires this to be **> tau** (strict) for
every transition, with mu = R9.4.1 6-mer mean levels in pA and tau = 1..8 pA. They
consider only the encoded strand; we additionally report the RC strand.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR = ROOT / "data" / "kmer_models"
CACHE_DIR = MODEL_DIR / "cache"
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b)] = _i


@dataclass(frozen=True)
class KmerModel:
    """Expected current level per k-mer, indexed by 2-bit code."""

    name: str
    k: int
    levels: np.ndarray  # float64, shape (4**k,)
    level_std: np.ndarray | None  # float64 or None
    units: str  # "std" or "pA"
    sha256: str = ""


def _parse_table(path: Path, k: int, has_std: bool) -> tuple[np.ndarray, np.ndarray | None]:
    n = 4**k
    levels = np.full(n, np.nan)
    std = np.full(n, np.nan) if has_std else None
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts or parts[0] == "kmer":
                continue
            c = int(kmer_code(parts[0]))
            levels[c] = float(parts[1])
            if std is not None:
                std[c] = float(parts[2])
    if np.isnan(levels).any() or (std is not None and np.isnan(std).any()):
        raise ValueError(f"{path} does not cover all {n} {k}-mers")
    return levels, std


_SPECS = {
    "r10.4.1_9mer": (9, False, "std"),
    "r9.4.1_6mer": (6, True, "pA"),
}


@lru_cache(maxsize=None)
def load_model(name: str) -> KmerModel:
    """Load a downloaded model, verifying its SHA-256 against PROVENANCE.json."""
    if name not in _SPECS:
        raise KeyError(f"unknown model {name!r}; choose from {sorted(_SPECS)}")
    prov_path = MODEL_DIR / "PROVENANCE.json"
    if not prov_path.exists():
        raise FileNotFoundError("k-mer tables missing: run scripts/fetch_kmer_models.py (never substitute a table)")
    entry = json.loads(prov_path.read_text())["files"][name]
    path = MODEL_DIR / entry["local_file"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != entry["sha256"]:
        raise ValueError(f"{path} sha256 {digest} != provenance {entry['sha256']}")
    k, has_std, units = _SPECS[name]
    levels, std = _parse_table(path, k, has_std)
    return KmerModel(name, k, levels, std, units, digest)


# --------------------------------------------------------------------------- codes


def kmer_code(kmer: str) -> int:
    """2-bit code of a k-mer string (MSB = first base)."""
    c = 0
    for b in kmer:
        c = (c << 2) | "ACGT".index(b)
    return c


def code_to_kmer(code: int, k: int) -> str:
    return "".join("ACGT"[(code >> (2 * (k - 1 - i))) & 3] for i in range(k))


def encode_kmers(seq: str, k: int) -> np.ndarray:
    """Codes of all ``len(seq)-k+1`` overlapping k-mers (int64), via the rolling update."""
    s = _LUT[np.frombuffer(seq.encode(), dtype=np.uint8)].astype(np.int64)
    if (s == 255).any():
        raise ValueError("sequence contains non-ACGT characters")
    n = len(s) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.int64)
    codes = np.zeros(n, dtype=np.int64)
    for i in range(k):  # vectorised equivalent of the rolling update
        codes = (codes << 2) | s[i : i + n]
    return codes


@lru_cache(maxsize=None)
def rc_code_table(k: int) -> np.ndarray:
    """``rc[c]`` = code of the reverse complement of k-mer ``c`` (complement = 3 - base)."""
    codes = np.arange(4**k, dtype=np.int64)
    out = np.zeros_like(codes)
    for i in range(k):
        base = (codes >> (2 * i)) & 3  # i-th base from the 3' end
        out = (out << 2) | (3 - base)
    return out


# --------------------------------------------------------------------------- confusability


def _delta_key(delta: float) -> str:
    return f"{delta:.6g}"


def confusability_table(model: KmerModel, delta: float | None = None, c: float | None = None, use_cache: bool = True) -> np.ndarray:
    """Per-k-mer confusability in [0, 1] (float32, shape (4**k,)).

    For ``units == "std"`` pass ``delta`` (default 0.05). For ``units == "pA"`` pass ``c``
    (default 1.0): delta(x, y) = c * sqrt((sd_x**2 + sd_y**2) / 2).
    """
    if model.units == "std":
        if c is not None:
            raise ValueError("std-unit model takes delta, not c")
        param = 0.05 if delta is None else float(delta)
        tag = f"delta{_delta_key(param)}"
    else:
        if delta is not None:
            raise ValueError("pA model takes c (multiple of pooled level std), not delta")
        param = 1.0 if c is None else float(c)
        tag = f"c{_delta_key(param)}"
    cache = CACHE_DIR / f"conf_{model.name}_{model.sha256[:12]}_{tag}.npy"
    if use_cache and cache.exists():
        return np.load(cache)
    k = model.k
    codes = np.arange(4**k, dtype=np.int64)
    count = np.zeros(4**k, dtype=np.int32)
    lv = model.levels
    for pos in range(k):
        shift = 2 * (k - 1 - pos)
        cur = (codes >> shift) & 3
        for off in (1, 2, 3):
            nb = codes ^ (((cur + off) & 3) ^ cur) << shift  # replace base at pos
            diff = np.abs(lv - lv[nb])
            if model.units == "std":
                thr = param
            else:
                sd = model.level_std
                thr = param * np.sqrt((sd**2 + sd[nb] ** 2) / 2.0)
            count += diff < thr
    table = (count / (3 * k)).astype(np.float32)
    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        np.save(cache, table)
    return table


def _summ(vals: np.ndarray) -> tuple[float, float]:
    return (float(vals.mean()), float(vals.max())) if vals.size else (float("nan"), float("nan"))


def oligo_confusability(seq: str, table: np.ndarray, k: int) -> dict[str, float]:
    """Mean and max confusability over the k-mers of ``seq`` and of its reverse complement.

    ``mean``/``max`` combine both strands (the mean over all k-mers of both strands; the max
    over both).
    """
    codes = encode_kmers(seq, k)
    fwd = table[codes]
    rc = table[rc_code_table(k)[codes]]
    fm, fx = _summ(fwd)
    rm, rx = _summ(rc)
    both = np.concatenate([fwd, rc])
    bm, bx = _summ(both)
    return {"mean": bm, "max": bx, "fwd_mean": fm, "fwd_max": fx, "rc_mean": rm, "rc_max": rx}


def adjacent_contrast(seq: str, model: KmerModel, tau: float | None = None) -> dict[str, float]:
    """Whritenour et al. adjacent-level contrast on both strands.

    ``tau`` defaults to 4 pA for pA models (their headline value) and 0.5 for std models
    (an arbitrary placeholder: there is no published tau in std units). ``frac_below_tau``
    is the fraction of transitions with ``|dl| <= tau``, i.e. violating their strict
    ``> tau`` constraint.
    """
    if tau is None:
        tau = 4.0 if model.units == "pA" else 0.5
    k = model.k
    codes = encode_kmers(seq, k)
    out: dict[str, float] = {"tau": float(tau)}
    for strand, cs in (("fwd", codes), ("rc", rc_code_table(k)[codes][::-1])):
        d = np.abs(np.diff(model.levels[cs]))
        out[f"{strand}_min"] = float(d.min()) if d.size else float("nan")
        out[f"{strand}_mean"] = float(d.mean()) if d.size else float("nan")
        out[f"{strand}_frac_below_tau"] = float((d <= tau).mean()) if d.size else float("nan")
    out["min"] = min(out["fwd_min"], out["rc_min"])
    out["mean"] = (out["fwd_mean"] + out["rc_mean"]) / 2
    out["frac_below_tau"] = (out["fwd_frac_below_tau"] + out["rc_frac_below_tau"]) / 2
    return out
