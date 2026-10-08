"""LT (Luby transform) fountain code: our own implementation.

Written from the descriptions in Luby (2002), "LT codes", FOCS 2002, and MacKay (2003/2005),
*Information Theory, Inference, and Learning Algorithms*, ch. 50. No code was taken from the
DNA Fountain codebase (Erlich & Zielinski 2017); only its published parameter *choices*
(c = 0.025, delta = 0.001) are reused as defaults, and they are cited as such.

Encoding
--------
The source is split into ``k`` equal-length segments. A droplet is identified by a 32-bit
seed. The seed seeds a PCG64 generator, which draws a degree ``d`` from the robust
soliton distribution and then ``d`` distinct neighbour indices (Floyd's algorithm, so it
is O(d) and not O(k)). The droplet payload is the XOR of those segments. The decoder
regenerates the neighbour set from the seed, so only ``(seed, payload)`` has to be stored.

Robust soliton (Luby 2002; MacKay eq. 50.4-50.6)
------------------------------------------------
``rho(1) = 1/k``, ``rho(d) = 1/(d(d-1))`` for ``d = 2..k``.
``S = c * ln(k/delta) * sqrt(k)``.
``tau(d) = S/(k d)`` for ``d < k/S``, ``tau(k/S) = S ln(S/delta)/k``, otherwise 0.
``mu = (rho + tau) / Z``.

Decoding
--------
1. Peeling (belief propagation): every droplet of residual degree 1 releases a segment,
   which is XORed out of every other droplet that contains it.
2. If peeling stalls, Gaussian elimination over GF(2) on the residual system (bit-packed
   coefficient rows, vectorised numpy XOR of byte rows) recovers whatever is
   determined.

Failure is always reported (``LTResult.ok = False``, ``n_missing > 0``). Corrupted droplets
are assumed to have been removed upstream (inner CRC); the decoder does not detect them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Iterable

import numpy as np

DEFAULT_C = 0.025  # Erlich & Zielinski (2017) parameter choice
DEFAULT_DELTA = 0.001  # Erlich & Zielinski (2017) parameter choice


class RobustSoliton:
    """Robust soliton degree distribution over ``1..k``."""

    def __init__(self, k: int, c: float = DEFAULT_C, delta: float = DEFAULT_DELTA) -> None:
        if k < 1:
            raise ValueError("k must be >= 1")
        self.k, self.c, self.delta = k, c, delta
        d = np.arange(1, k + 1, dtype=np.float64)
        rho = np.zeros(k)
        rho[0] = 1.0 / k
        if k > 1:
            rho[1:] = 1.0 / (d[1:] * (d[1:] - 1.0))
        tau = np.zeros(k)
        S = c * math.log(k / delta) * math.sqrt(k)
        if S > 0:
            spike = int(round(k / S))
            spike = min(max(spike, 1), k)
            for i in range(1, spike):
                tau[i - 1] = S / (k * i)
            tau[spike - 1] += S * math.log(S / delta) / k if S > delta else 0.0
        mu = rho + np.maximum(tau, 0.0)
        self.Z = float(mu.sum())
        self.pmf = mu / self.Z
        self.cdf = np.cumsum(self.pmf)
        self.cdf[-1] = 1.0
        self.S = S

    def sample(self, rng: np.random.Generator) -> int:
        return int(np.searchsorted(self.cdf, rng.random(), side="right")) + 1

    @property
    def mean_degree(self) -> float:
        return float((np.arange(1, self.k + 1) * self.pmf).sum())


@dataclass(frozen=True)
class Droplet:
    seed: int
    data: bytes


def droplet_neighbors(seed: int, k: int, dist: RobustSoliton) -> list[int]:
    """Deterministic, sorted neighbour set of the droplet with 32-bit ``seed``."""
    if not 0 <= seed < 2**32:
        raise ValueError("seed must be a 32-bit unsigned integer")
    if dist.k != k:
        raise ValueError("distribution k mismatch")
    rng = np.random.Generator(np.random.PCG64(seed))
    d = min(dist.sample(rng), k)
    chosen: set[int] = set()
    for j in range(k - d, k):  # Floyd's algorithm: d distinct values from range(k)
        t = int(rng.integers(0, j + 1))
        chosen.add(j if t in chosen else t)
    return sorted(chosen)


def _xor_segments(segments: list[bytes], idx: list[int]) -> bytes:
    acc = np.frombuffer(segments[idx[0]], dtype=np.uint8).copy()
    for i in idx[1:]:
        acc ^= np.frombuffer(segments[i], dtype=np.uint8)
    return acc.tobytes()


def encode_droplet(segments: list[bytes], seed: int, dist: RobustSoliton) -> bytes:
    """XOR of the droplet's neighbour segments (all segments must have equal length)."""
    return _xor_segments(segments, droplet_neighbors(seed, len(segments), dist))


def generate_screened(
    segments: list[bytes],
    seed_iter: Iterable[int],
    accept: Callable[[int, bytes], bool],
    n_needed: int,
    dist: RobustSoliton | None = None,
    stats: dict | None = None,
) -> list[tuple[int, bytes]]:
    """Draw droplets from ``seed_iter`` and keep those for which ``accept(seed, data)`` holds.

    Stops when ``n_needed`` droplets are accepted. Raises ``RuntimeError`` if the seed
    iterator runs out first (never returns fewer droplets silently). Counts go into
    ``stats`` (``drawn``, ``accepted``, ``rejected``) when it is provided.
    """
    dist = dist or RobustSoliton(len(segments))
    out: list[tuple[int, bytes]] = []
    drawn = rejected = 0
    for seed in seed_iter:
        if len(out) >= n_needed:
            break
        drawn += 1
        data = encode_droplet(segments, seed, dist)
        if accept(seed, data):
            out.append((seed, data))
        else:
            rejected += 1
    if stats is not None:
        stats.update(drawn=drawn, accepted=len(out), rejected=rejected)
    if len(out) < n_needed:
        raise RuntimeError(f"seed iterator exhausted: {len(out)}/{n_needed} droplets accepted")
    return out


@dataclass
class LTResult:
    ok: bool
    segments: list[bytes] | None
    n_missing: int
    n_peeled: int
    n_ge: int
    n_droplets: int
    n_duplicates: int
    info: dict = field(default_factory=dict)


class LTDecoder:
    """Collects droplets, then decodes by peeling with a GF(2) Gaussian-elimination fallback."""

    def __init__(self, k: int, seg_len: int, dist: RobustSoliton | None = None, max_ge_unknowns: int = 20000) -> None:
        self.k, self.seg_len = k, seg_len
        self.dist = dist or RobustSoliton(k)
        self.max_ge_unknowns = max_ge_unknowns
        self._seeds: set[int] = set()
        self._droplets: list[tuple[list[int], np.ndarray]] = []
        self.n_duplicates = 0

    def add(self, seed: int, data: bytes) -> None:
        if len(data) != self.seg_len:
            raise ValueError("droplet length mismatch")
        if seed in self._seeds:
            self.n_duplicates += 1
            return
        self._seeds.add(seed)
        self._droplets.append((droplet_neighbors(seed, self.k, self.dist), np.frombuffer(data, dtype=np.uint8).copy()))

    def decode(self) -> LTResult:
        k = self.k
        nbrs = [set(n) for n, _ in self._droplets]
        vals = [v.copy() for _, v in self._droplets]
        seg_of: list[list[int]] = [[] for _ in range(k)]
        for di, ns in enumerate(nbrs):
            for j in ns:
                seg_of[j].append(di)
        solved: list[np.ndarray | None] = [None] * k
        ripple = [di for di, ns in enumerate(nbrs) if len(ns) == 1]
        n_peeled = 0
        while ripple:
            di = ripple.pop()
            if len(nbrs[di]) != 1:
                continue
            j = next(iter(nbrs[di]))
            if solved[j] is not None:
                nbrs[di].clear()
                continue
            solved[j] = vals[di]
            n_peeled += 1
            for dj in seg_of[j]:
                if j in nbrs[dj]:
                    nbrs[dj].discard(j)
                    if dj != di:
                        vals[dj] = vals[dj] ^ solved[j]
                    if len(nbrs[dj]) == 1:
                        ripple.append(dj)
        unknown = [j for j in range(k) if solved[j] is None]
        n_ge = 0
        info: dict = {}
        if unknown and len(unknown) <= self.max_ge_unknowns:
            n_ge = self._gaussian(unknown, nbrs, vals, solved, info)
        missing = sum(s is None for s in solved)
        segs = [s.tobytes() for s in solved] if missing == 0 else None
        return LTResult(missing == 0, segs, missing, n_peeled, n_ge, len(self._droplets), self.n_duplicates, info)

    @staticmethod
    def _gaussian(unknown: list[int], nbrs: list[set[int]], vals: list[np.ndarray], solved: list, info: dict) -> int:
        """GF(2) elimination on residual droplets restricted to ``unknown`` columns."""
        col = {j: c for c, j in enumerate(unknown)}
        u = len(unknown)
        rows = [di for di, ns in enumerate(nbrs) if ns]
        info["ge_unknowns"], info["ge_rows"] = u, len(rows)
        if not rows:
            return 0
        nbytes = (u + 7) // 8
        A = np.zeros((len(rows), nbytes), dtype=np.uint8)
        B = np.stack([vals[di] for di in rows])
        for r, di in enumerate(rows):
            for j in nbrs[di]:
                c = col[j]  # every residual neighbour is an unknown by construction of peeling
                A[r, c >> 3] |= np.uint8(0x80 >> (c & 7))
        m = len(rows)
        pivot_row_of_col = [-1] * u
        r = 0
        for c in range(u):
            if r >= m:
                break
            byte, mask = c >> 3, np.uint8(0x80 >> (c & 7))
            colbits = (A[r:, byte] & mask) != 0
            hits = np.nonzero(colbits)[0]
            if hits.size == 0:
                continue
            p = r + int(hits[0])
            if p != r:
                A[[r, p]] = A[[p, r]]
                B[[r, p]] = B[[p, r]]
            others = np.nonzero((A[:, byte] & mask) != 0)[0]
            others = others[others != r]
            if others.size:
                A[others] ^= A[r]
                B[others] ^= B[r]
            pivot_row_of_col[c] = r
            r += 1
        n = 0
        for c, pr in enumerate(pivot_row_of_col):
            if pr < 0:
                continue
            # solved iff the pivot row has no other unknown left (fully reduced form)
            row = A[pr].copy()
            row[c >> 3] &= np.uint8(~(0x80 >> (c & 7)) & 0xFF)
            if not row.any():
                solved[unknown[c]] = B[pr].copy()
                n += 1
        return n
