"""Low-level sequence utilities shared by all codecs.

Conventions
-----------
* Sequences are upper-case ``str`` over ``ACGT``.
* Rotating ternary (Goldman et al. 2013): given the previous base ``p`` and a trit
  ``t`` in {0,1,2}, the next base is ``BASES[(BASES.index(p) + 1 + t) % 4]``.
  The next base therefore never equals ``p``, so runs have length 1.
* Byte <-> trit conversion is a per-frame big-integer conversion with a fixed
  number of trits (DESIGN.md F14).
"""

from __future__ import annotations

import math

import numpy as np

BASES = "ACGT"
_BASE_IDX = {b: i for i, b in enumerate(BASES)}
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    """Reverse complement of a DNA string."""
    return seq.translate(_COMP)[::-1]


# --------------------------------------------------------------------------- 2-bit


def bytes_to_bases_2bit(data: bytes) -> str:
    """Map bytes to bases, 2 bits per base, MSB first (00->A, 01->C, 10->G, 11->T)."""
    out = []
    for byte in data:
        for shift in (6, 4, 2, 0):
            out.append(BASES[(byte >> shift) & 3])
    return "".join(out)


def bases_to_bytes_2bit(seq: str) -> bytes:
    """Inverse of :func:`bytes_to_bases_2bit`. ``len(seq)`` must be a multiple of 4."""
    if len(seq) % 4:
        raise ValueError("2-bit sequence length must be a multiple of 4")
    out = bytearray()
    for i in range(0, len(seq), 4):
        v = 0
        for b in seq[i : i + 4]:
            v = (v << 2) | _BASE_IDX[b]
        out.append(v)
    return bytes(out)


# --------------------------------------------------------------------------- ternary


def max_bytes_for_trits(n_trits: int) -> int:
    """Largest byte count ``nb`` with ``256**nb <= 3**n_trits`` (exact integer check)."""
    if n_trits <= 0:
        return 0
    nb = int(n_trits * math.log2(3) // 8)
    while 256 ** (nb + 1) <= 3**n_trits:
        nb += 1
    while nb > 0 and 256**nb > 3**n_trits:
        nb -= 1
    return nb


def min_trits_for_bytes(n_bytes: int) -> int:
    """Smallest trit count able to carry ``n_bytes`` bytes."""
    t = 0
    while 3**t < 256**n_bytes:
        t += 1
    return t


def bytes_to_trits(data: bytes, n_trits: int) -> list[int]:
    """Big-endian big-integer conversion of ``data`` to exactly ``n_trits`` trits."""
    if 256 ** len(data) > 3**n_trits:
        raise ValueError(f"{len(data)} bytes do not fit in {n_trits} trits")
    v = int.from_bytes(data, "big")
    trits = [0] * n_trits
    for i in range(n_trits - 1, -1, -1):
        v, trits[i] = divmod(v, 3)
    return trits


def trits_to_bytes(trits: list[int], n_bytes: int) -> bytes | None:
    """Inverse of :func:`bytes_to_trits`. Returns ``None`` if the value overflows ``n_bytes``.

    An overflow can only arise from a corrupted sequence, so it is a free error check.
    """
    v = 0
    for t in trits:
        v = v * 3 + t
    if v >= 256**n_bytes:
        return None
    return v.to_bytes(n_bytes, "big")


def int_to_trits(value: int, n_trits: int) -> list[int]:
    """Fixed-width base-3 representation of a non-negative integer."""
    if value < 0 or value >= 3**n_trits:
        raise ValueError(f"{value} does not fit in {n_trits} trits")
    trits = [0] * n_trits
    for i in range(n_trits - 1, -1, -1):
        value, trits[i] = divmod(value, 3)
    return trits


def trits_to_int(trits: list[int]) -> int:
    v = 0
    for t in trits:
        v = v * 3 + t
    return v


def rotate_encode(trits: list[int], prev: str) -> str:
    """Rotating-ternary encoding of ``trits`` starting after base ``prev``."""
    out = []
    p = _BASE_IDX[prev]
    for t in trits:
        p = (p + 1 + t) % 4
        out.append(BASES[p])
    return "".join(out)


def rotate_decode(seq: str, prev: str) -> list[int] | None:
    """Inverse of :func:`rotate_encode`. Returns ``None`` if a base repeats its predecessor
    (impossible in a valid codeword, so it signals an error)."""
    out = []
    p = _BASE_IDX[prev]
    for b in seq:
        q = _BASE_IDX[b]
        t = (q - p - 1) % 4
        if t == 3:
            return None
        out.append(t)
        p = q
    return out


# --------------------------------------------------------------------------- PRNG


def keystream(n: int, *key: int) -> bytes:
    """Deterministic pseudo-random bytes keyed by a tuple of non-negative ints.

    Uses numpy's PCG64 seeded through ``SeedSequence``, which is stable across platforms
    for a fixed numpy version (pinned in pyproject.toml).
    """
    if n == 0:
        return b""
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(list(key))))
    return rng.bytes(n)


def xor_bytes(a: bytes, b: bytes) -> bytes:
    """Bytewise XOR of equal-length byte strings."""
    if len(a) != len(b):
        raise ValueError("length mismatch")
    return (np.frombuffer(a, np.uint8) ^ np.frombuffer(b, np.uint8)).tobytes()


# --------------------------------------------------------------------------- bit frames


def max_bits_for_trits(n_trits: int) -> int:
    """Largest ``b`` with ``2**b <= 3**n_trits``."""
    return (3**n_trits).bit_length() - 1 if n_trits > 0 else 0


def keystream_int(nbits: int, *key: int) -> int:
    """``nbits`` pseudo-random bits as an int (see :func:`keystream`)."""
    if nbits <= 0:
        return 0
    nb = -(-nbits // 8)
    return int.from_bytes(keystream(nb, *key), "big") >> (8 * nb - nbits)


def int_to_base4(value: int, n: int) -> str:
    """Fixed-width base-4 digits of ``value`` (MSB first) as bases A/C/G/T = 0/1/2/3."""
    if value < 0 or value >> (2 * n):
        raise ValueError("value does not fit")
    return "".join(BASES[(value >> (2 * (n - 1 - i))) & 3] for i in range(n))


def base4_to_int(seq: str) -> int:
    v = 0
    for b in seq:
        v = (v << 2) | _BASE_IDX[b]
    return v
