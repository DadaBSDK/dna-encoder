"""Noise channel simulator: coverage, dropout, strand orientation and IDS errors.

Model (exact, see docs/notes/channel.md)
----------------------------------------
1. **Copy weights.** Each oligo ``i`` gets a weight ``w_i``:

   * ``coverage="poisson"``: ``w_i = 1`` (reads are Poisson-distributed around the mean);
   * ``coverage="lognormal"``: ``w_i ~ LogNormal(-sigma**2/2, sigma)`` (mean 1);
   * ``coverage="model"``: ``w_i = dnastore.coverage.expected_copies(features, preset)``
     (synthesis + PCR yield, sources in docs/notes/pcr_bias.md). ``features`` must hold
     ``"gc"`` and may hold ``"accessibility"``. ``zero_accessibility=True`` sets the
     ASSUMED accessibility penalty to 0. Every coverage-model result must also be
     reported with this switch on (user requirement).

2. **Dropout.** Each oligo is lost independently with probability ``dropout`` (``w_i = 0``).
3. **Read sampling.** The total read count is fixed at ``N = round(mean_coverage * n_oligos)``
   (a sequencing run of fixed depth) and allocated by ``Multinomial(N, w / sum(w))``. Dropout
   therefore redistributes depth to the surviving oligos rather than lowering the total.
4. **Orientation.** Each read is the reverse complement of its oligo with probability
   ``rc_frac``. Errors are applied **in read orientation**, i.e. to the strand as sequenced.
5. **IDS errors**, independently per template position ``j`` of the read-orientation
   strand. With ``u ~ U(0,1)``: delete if ``u < p_del,j``; else substitute (uniformly by one of
   the 3 other bases) if ``u < p_del,j + p_sub,j``. Independently, one uniformly random base
   is inserted *before* position ``j`` with probability ``p_ins,j``.

6. **Confusability-weighted mode** (``conf_weighted=True``). This is a **PROXY** for
   nanopore behaviour, **NOT a validated error model**. It is also **CIRCULAR** with respect
   to the steering encoder's confusability objective (it rewards exactly what the encoder
   optimises), so it is reported only as a sanity check. Rates become ``p_x,j = p_x * m_j``
   with ``m_j = cbar_j / mu``, where ``cbar_j`` is the mean confusability of the k-mers
   (fully inside the strand) that cover position ``j``, and ``mu`` is the mean of the
   per-k-mer table. For a uniform random sequence every covering k-mer is uniform over all
   k-mers, so ``E[m_j] = 1`` exactly, and the expected per-base rate equals the base rate.
   Rates are clipped to ``p_del + p_sub <= 0.9`` and ``p_ins <= 0.9``. Multipliers are
   computed on the read-orientation strand (for an RC read, on ``revcomp(oligo)``).
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
import math
from typing import Any

import numpy as np

_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


@dataclass
class ChannelParams:
    p_sub: float = 0.0
    p_ins: float = 0.0
    p_del: float = 0.0
    dropout: float = 0.0
    mean_coverage: float = 10.0
    coverage: str = "poisson"  # "poisson" | "lognormal" | "model"
    lognormal_sigma: float = 0.3
    coverage_preset: str = "twist_pcr_measured"  # key of dnastore.coverage.SENSITIVITY_GRID
    zero_accessibility: bool = False  # set the ASSUMED accessibility knob to 0
    rc_frac: float = 0.5
    conf_weighted: bool = False
    conf_model: str = "r10.4.1_9mer"
    conf_delta: float = 0.05
    seed: int = 0

    def __post_init__(self) -> None:
        for name in ("p_sub", "p_ins", "p_del", "dropout", "rc_frac"):
            v = getattr(self, name)
            if not math.isfinite(v) or not 0 <= v <= 1:
                raise ValueError(f"{name} must be in [0, 1]")
        if self.p_sub + self.p_del > 1:
            raise ValueError("p_sub + p_del must not exceed 1")
        for name in ("mean_coverage", "lognormal_sigma"):
            v = getattr(self, name)
            if not math.isfinite(v) or v < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.coverage not in ("poisson", "lognormal", "model"):
            raise ValueError(f"unknown coverage mode {self.coverage!r}")

    @classmethod
    def from_dict(cls, d: dict[str, Any] | None) -> ChannelParams:
        d = dict(d or {})
        unknown = set(d) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"unknown ChannelParams keys: {sorted(unknown)}")
        return cls(**d)


@dataclass
class Read:
    seq: str
    oligo: int  # index into the input oligo list
    reverse: bool
    n_sub: int = 0
    n_ins: int = 0
    n_del: int = 0


@dataclass
class SimulationTruth:
    """Per-oligo bookkeeping of a simulation."""

    weights: np.ndarray
    dropped: np.ndarray
    read_counts: np.ndarray
    extra: dict[str, Any] = field(default_factory=dict)


def _encode(seqs: list[str]) -> np.ndarray:
    L = len(seqs[0])
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), L)
    return _LUT[arr]


def _revcomp_codes(a: np.ndarray) -> np.ndarray:
    return (3 - a)[..., ::-1]


def copy_weights(oligos: list[str], params: ChannelParams, rng: np.random.Generator,
                 features: dict[str, Any] | None = None) -> np.ndarray:
    """Unnormalised per-oligo copy weights (before dropout)."""
    n = len(oligos)
    if params.coverage == "poisson":
        return np.ones(n)
    if params.coverage == "lognormal":
        s = params.lognormal_sigma
        return rng.lognormal(-0.5 * s * s, s, n)
    if params.coverage == "model":
        from . import coverage as cov

        model = cov.SENSITIVITY_GRID[params.coverage_preset]
        if params.zero_accessibility:
            model = cov.with_overrides(model, acc_penalty=0.0)
        feats = dict(features or {})
        if "gc" not in feats:
            feats["gc"] = np.array([(s.count("G") + s.count("C")) / len(s) for s in oligos])
        if params.zero_accessibility:
            feats.pop("accessibility", None)
        return np.asarray(cov.expected_copies(feats, model, rng), dtype=float)
    raise ValueError(f"unknown coverage mode {params.coverage!r}")


_CONF_CACHE: dict[tuple, tuple[np.ndarray, int, float]] = {}


def _conf(params: ChannelParams) -> tuple[np.ndarray, int, float]:
    key = (params.conf_model, params.conf_delta)
    if key not in _CONF_CACHE:
        from . import kmer

        m = kmer.load_model(params.conf_model)
        t = kmer.confusability_table(m, delta=params.conf_delta) if m.units == "std" else kmer.confusability_table(m, c=params.conf_delta)
        t = np.asarray(t, dtype=np.float64)
        _CONF_CACHE[key] = (t, m.k, float(t.mean()))
    return _CONF_CACHE[key]


def conf_multipliers(seq: str, params: ChannelParams) -> np.ndarray:
    """Per-position error multipliers ``m_j = cbar_j / mu`` (see module docstring)."""
    from . import kmer

    table, k, mu = _conf(params)
    L = len(seq)
    if L < k:
        return np.ones(L)
    c = table[kmer.encode_kmers(seq, k)]  # conf of k-mer starting at s, s = 0..L-k
    cs = np.concatenate([[0.0], np.cumsum(c)])
    j = np.arange(L)
    lo = np.clip(j - k + 1, 0, L - k)  # first covering k-mer start
    hi = np.clip(j, 0, L - k)  # last covering k-mer start
    return (cs[hi + 1] - cs[lo]) / (hi - lo + 1) / mu


def _mutate_batch(tmpl: np.ndarray, rates: tuple[np.ndarray, np.ndarray, np.ndarray],
                  rng: np.random.Generator) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Apply IDS errors to a (n, L) batch of base codes. ``rates`` are (n, L) arrays or scalars."""
    n, L = tmpl.shape
    p_sub, p_ins, p_del = rates
    u = rng.random((n, L))
    is_del = u < p_del
    is_sub = (~is_del) & (u < p_del + p_sub)
    is_ins = rng.random((n, L)) < p_ins
    body = tmpl.copy()
    body[is_sub] = (body[is_sub] + rng.integers(1, 4, size=int(is_sub.sum()), dtype=np.uint8)) % 4
    body[is_del] = 255
    ins = np.full((n, L), 255, dtype=np.uint8)
    ins[is_ins] = rng.integers(0, 4, size=int(is_ins.sum()), dtype=np.uint8)
    tok = np.stack([ins, body], axis=2).reshape(n, 2 * L)
    keep = tok != 255
    lens = keep.sum(axis=1)
    flat = _BASES[tok[keep]].tobytes().decode()
    ends = np.cumsum(lens)
    starts = ends - lens
    seqs = [flat[s:e] for s, e in zip(starts, ends)]
    return seqs, is_sub.sum(1), is_ins.sum(1), is_del.sum(1)


def simulate(oligos: list[str], params: ChannelParams, rng: np.random.Generator | None = None,
             features: dict[str, Any] | None = None, batch: int = 50_000) -> tuple[list[Read], SimulationTruth]:
    """Simulate a sequencing run over ``oligos`` (all of equal length). Returns (reads, truth)."""
    if rng is None:
        rng = np.random.default_rng(params.seed)
    n = len(oligos)
    if n == 0:
        return [], SimulationTruth(np.zeros(0), np.zeros(0, bool), np.zeros(0, int))
    if len({len(s) for s in oligos}) != 1:
        raise ValueError("simulate() expects equal-length oligos")
    if batch < 1 or any(not s or set(s) - set("ACGT") for s in oligos):
        raise ValueError("simulate() requires nonempty ACGT oligos and a positive batch size")
    w = copy_weights(oligos, params, rng, features)
    dropped = rng.random(n) < params.dropout
    w = np.where(dropped, 0.0, w)
    total = int(round(params.mean_coverage * n))
    counts = rng.multinomial(total, w / w.sum()) if w.sum() > 0 and total > 0 else np.zeros(n, dtype=int)

    codes = _encode(oligos)
    src = np.repeat(np.arange(n), counts)
    rng.shuffle(src)
    rev = rng.random(src.size) < params.rc_frac

    mults = None
    if params.conf_weighted:
        mults = np.stack([np.stack([conf_multipliers(s, params) for s in oligos]),
                          np.stack([conf_multipliers(_rc(s), params) for s in oligos])])  # (2, n, L)

    reads: list[Read] = []
    for b0 in range(0, src.size, batch):
        s_idx, r_flag = src[b0 : b0 + batch], rev[b0 : b0 + batch]
        tmpl = codes[s_idx]
        if r_flag.any():
            tmpl[r_flag] = _revcomp_codes(tmpl[r_flag])
        if mults is None:
            rates = (params.p_sub, params.p_ins, params.p_del)
        else:
            m = mults[r_flag.astype(int), s_idx]
            pd, ps = params.p_del * m, params.p_sub * m
            scale = np.minimum(1.0, 0.9 / np.maximum(pd + ps, 1e-300))
            rates = (ps * scale, np.minimum(params.p_ins * m, 0.9), pd * scale)
        seqs, ns, ni, nd = _mutate_batch(tmpl, rates, rng)
        reads.extend(Read(sq, int(o), bool(r), int(a), int(b), int(c))
                     for sq, o, r, a, b, c in zip(seqs, s_idx, r_flag, ns, ni, nd))
    return reads, SimulationTruth(w, dropped, counts)


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
