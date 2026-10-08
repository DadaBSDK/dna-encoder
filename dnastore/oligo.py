"""Oligo assembly: primers, the shared index field, and header oligos.

Physical layout (DESIGN.md §2.1, Part 0.8)::

    [F primer][index: INDEX_TRITS nt][seed: SEED_TRITS nt][payload region][revcomp(R primer)]

* The **index field** is codec-independent: rotating ternary chained from the last base of
  F, so it has no homopolymers and the junction with F has none either. Every read can be
  grouped and its header/data role identified before the codec is known.
* The **seed field** (SEED_TRITS rotating-ternary nt, also codec-independent) selects one of
  ``3**SEED_TRITS`` re-encodings of the same frame. The encoder rerolls the seed when an
  oligo violates a framework hard constraint (3'-anchored primer match, for ALL arms) or a
  codec's own hard constraints. How a codec uses the seed is codec-specific: e.g. seed 0 =
  identity and seed > 0 = keystream whitening for the baselines. The seed has to be readable
  before any keystream is known, so it is not whitened (DESIGN.md F7).
* Indices ``0 .. HEADER_INDEX_SPACE-1`` are reserved for header fragments
  (index = fragment * copies + copy). Data and parity oligos start at ``HEADER_INDEX_SPACE``.
* **Header oligos** carry ``[frag_id<<4 | n_frags][chunk]`` plus CRC-32, XOR-whitened with a
  keystream keyed by the fixed ``HEADER_KEY`` and the oligo index. The R copies of a
  fragment therefore have completely different sequences and fail independently
  (DESIGN.md F6). Their payload region is rotating ternary, not steered.
* There is no orientation marker: orientation comes from the asymmetric primers (F8).
"""

from __future__ import annotations

from dataclasses import dataclass

from .dnautil import (
    bytes_to_trits,
    int_to_trits,
    keystream,
    max_bytes_for_trits,
    rotate_decode,
    rotate_encode,
    trits_to_bytes,
    trits_to_int,
    xor_bytes,
)
from .ecc import CRC_BYTES, add_crc, check_crc

INDEX_TRITS = 15  # 3**15 = 14,348,907 indices
SEED_TRITS = 3  # default seed field: 27 seeds. Data oligos may use 4 (81 seeds), stored in the header.
HEADER_SEED_TRITS = 3  # header oligos always use 3 so they can be parsed before the header is known
N_SEEDS = 3**SEED_TRITS
PREFIX_NT = INDEX_TRITS + SEED_TRITS
HEADER_INDEX_SPACE = 64
HEADER_KEY = 0x48445231  # "HDR1"; fixed, codec-independent whitening key for header oligos
MAX_HEADER_FRAGMENTS = 15


@dataclass(frozen=True)
class OligoGeometry:
    """Lengths derived from the total oligo length and the primers."""

    oligo_len: int
    f_len: int
    r_len: int
    seed_trits: int = SEED_TRITS

    @property
    def body_nt(self) -> int:
        return self.oligo_len - self.f_len - self.r_len

    @property
    def prefix_nt(self) -> int:
        return INDEX_TRITS + self.seed_trits

    @property
    def payload_nt(self) -> int:
        return self.body_nt - self.prefix_nt

    @property
    def header_payload_nt(self) -> int:
        return self.body_nt - INDEX_TRITS - HEADER_SEED_TRITS

    @property
    def n_seeds(self) -> int:
        return 3**self.seed_trits

    def validate(self) -> None:
        if not 1 <= self.seed_trits <= 8:
            raise ValueError("seed_trits must be in [1, 8]")
        if self.payload_nt < 16:
            raise ValueError(f"oligo_len {self.oligo_len} leaves only {self.payload_nt} payload nt")


def encode_index(index: int, prev: str) -> str:
    """Index field: ``INDEX_TRITS`` rotating-ternary bases following base ``prev``."""
    return rotate_encode(int_to_trits(index, INDEX_TRITS), prev)


def decode_index(seq: str, prev: str) -> int | None:
    """Decode an index field (``None`` if it contains an impossible repeat)."""
    if len(seq) != INDEX_TRITS:
        return None
    trits = rotate_decode(seq, prev)
    return None if trits is None else trits_to_int(trits)


def encode_seed(seed: int, prev: str, seed_trits: int = SEED_TRITS) -> str:
    """Seed field: ``seed_trits`` rotating-ternary bases following base ``prev``."""
    return rotate_encode(int_to_trits(seed, seed_trits), prev)


def decode_seed(seq: str, prev: str) -> int | None:
    trits = rotate_decode(seq, prev)
    return None if trits is None else trits_to_int(trits)


def encode_prefix(index: int, seed: int, prev: str, seed_trits: int = SEED_TRITS) -> str:
    """Index field followed by seed field."""
    idx = encode_index(index, prev)
    return idx + encode_seed(seed, idx[-1], seed_trits)


def decode_prefix(body: str, prev: str, seed_trits: int = SEED_TRITS) -> tuple[int, int] | None:
    """Return ``(index, seed)`` from the start of a body, or ``None``."""
    idx = decode_index(body[:INDEX_TRITS], prev)
    if idx is None:
        return None
    seed = decode_seed(body[INDEX_TRITS : INDEX_TRITS + seed_trits], body[INDEX_TRITS - 1])
    return None if seed is None else (idx, seed)


def assemble(forward: str, body: str, tail: str) -> str:
    """Full oligo ``F + body + revcomp(R)`` (``tail`` is already reverse-complemented)."""
    return forward + body + tail


# --------------------------------------------------------------------------- header oligos


def header_chunk_capacity(payload_nt: int) -> int:
    """Header bytes per fragment: ternary capacity minus CRC and the fragment byte."""
    return max_bytes_for_trits(payload_nt) - CRC_BYTES - 1


def header_fragments(header: bytes, payload_nt: int) -> list[bytes]:
    """Split packed header bytes into fixed-size fragment frames (before CRC)."""
    cap = header_chunk_capacity(payload_nt)
    if cap <= 0:
        raise ValueError("payload region too short for header fragments")
    chunks = [header[i : i + cap] for i in range(0, len(header), cap)] or [b""]
    n = len(chunks)
    if n > MAX_HEADER_FRAGMENTS:
        raise ValueError(f"header needs {n} fragments (> {MAX_HEADER_FRAGMENTS})")
    return [bytes([(i << 4) | n]) + c + bytes(cap - len(c)) for i, c in enumerate(chunks)]


def header_index(fragment: int, copy: int, copies: int) -> int:
    idx = fragment * copies + copy
    if idx >= HEADER_INDEX_SPACE:
        raise ValueError("too many header oligos for the reserved index space")
    return idx


def encode_header_payload(frag: bytes, index: int, seed: int, prev: str, payload_nt: int) -> str:
    """Rotating-ternary payload of one header oligo (whitened by (index, seed), CRC-protected)."""
    frame = add_crc(index, frag)
    frame = xor_bytes(frame, keystream(len(frame), HEADER_KEY, index, seed))
    return rotate_encode(bytes_to_trits(frame, payload_nt), prev)


def decode_header_payload(seq: str, index: int, seed: int, prev: str) -> tuple[int, int, bytes] | None:
    """Return ``(frag_id, n_frags, chunk)`` or ``None`` if invalid."""
    trits = rotate_decode(seq, prev)
    if trits is None:
        return None
    nb = max_bytes_for_trits(len(seq))
    frame = trits_to_bytes(trits, nb)
    if frame is None:
        return None
    frame = xor_bytes(frame, keystream(len(frame), HEADER_KEY, index, seed))
    payload = check_crc(index, frame)
    if payload is None or not payload:
        return None
    return payload[0] >> 4, payload[0] & 0x0F, payload[1:]
