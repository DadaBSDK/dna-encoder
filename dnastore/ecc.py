"""Error correction: outer Reed-Solomon across oligos, inner CRC-32 per oligo.

Outer code
----------
A systematic Reed-Solomon code over GF(256) (primitive polynomial 0x11D, alpha = 2,
first consecutive root alpha**0), defined through its parity-check matrix
``H[i, j] = alpha**(i*(n-1-j))`` for ``i < p`` and ``j < n = k + p <= 255``. Symbol ``j`` is
the coefficient of ``x**(n-1-j)``, the same convention as ``reedsolo`` (message first,
parity last, ``fcr=0``), so codewords are bit-identical to
``reedsolo.RSCodec(p, fcr=0, prim=0x11d, generator=2)`` (tested). Any ``p`` columns of ``H``
form a Vandermonde matrix on distinct nodes, so the code is MDS: any ``p`` erased symbols
are recoverable.

The data are arranged as a matrix of *rows* (``row_bytes``; see RSLayout for the bit packing
bytes) and each *column* is one RS codeword running across the ``k + p`` oligos of a
block. An oligo dropout erases the same row in every column of its block, so erasure
decoding is a single ``e x e`` solve per block applied to all columns at once. This is
the reason for a custom numpy implementation instead of per-codeword ``reedsolo``
(measured at ~6 ms encode / ~18 ms decode per codeword in pure Python, which is too slow
for the robustness sweeps).

Undetected inner errors (a corrupted oligo that passes the CRC, probability about 2**-32
per tested candidate) are handled by a consistency check on the spare syndromes plus a
search for a single bad row. More than one undetected bad row in a block is not
corrected, and the block is then reported as failed or inconsistent, never silently
accepted.

Inner check
-----------
CRC-32 (IEEE 802.3, ``zlib.crc32``) over ``index (4 bytes BE) || payload``. Binding the
index into the CRC means a read whose index was corrupted into another valid index fails
the check. It was CRC-16 until format version 3. The decoder's CRC-guided repair tests up
to ~10^4 candidate bodies per failed oligo, which would make CRC-16 false accepts
non-negligible (~0.15 per oligo) but leaves CRC-32 at ~2e-6.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------- GF(256)

_PRIM = 0x11D
GF_EXP = np.zeros(512, dtype=np.uint8)
GF_LOG = np.zeros(256, dtype=np.int32)
_x = 1
for _i in range(255):
    GF_EXP[_i] = _x
    GF_LOG[_x] = _i
    _x <<= 1
    if _x & 0x100:
        _x ^= _PRIM
GF_EXP[255:510] = GF_EXP[0:255]

_a = np.arange(256)
GF_MUL = np.zeros((256, 256), dtype=np.uint8)
_nz = _a[1:]
GF_MUL[1:, 1:] = GF_EXP[(GF_LOG[_nz][:, None] + GF_LOG[_nz][None, :]) % 255]
GF_INV = np.zeros(256, dtype=np.uint8)
GF_INV[1:] = GF_EXP[(255 - GF_LOG[_nz]) % 255]


def gf_matmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Matrix product over GF(256): ``a`` is (r, t), ``b`` is (t, c)."""
    a = np.asarray(a, dtype=np.uint8)
    b = np.asarray(b, dtype=np.uint8)
    out = np.zeros((a.shape[0], b.shape[1]), dtype=np.uint8)
    for t in range(a.shape[1]):
        out ^= GF_MUL[a[:, t][:, None], b[t][None, :]]
    return out


def gf_inv_matrix(m: np.ndarray) -> np.ndarray:
    """Inverse of a square matrix over GF(256) by Gauss-Jordan elimination."""
    n = m.shape[0]
    aug = np.concatenate([np.asarray(m, dtype=np.uint8), np.eye(n, dtype=np.uint8)], axis=1)
    for col in range(n):
        piv = next((r for r in range(col, n) if aug[r, col] != 0), None)
        if piv is None:
            raise np.linalg.LinAlgError("singular matrix over GF(256)")
        if piv != col:
            aug[[col, piv]] = aug[[piv, col]]
        aug[col] = GF_MUL[GF_INV[aug[col, col]], aug[col]]
        for r in range(n):
            if r != col and aug[r, col]:
                aug[r] ^= GF_MUL[aug[r, col], aug[col]]
    return aug[:, n:]


def parity_check_matrix(n: int, p: int) -> np.ndarray:
    """``H[i, j] = alpha**(i*(n-1-j))`` for ``i < p``, ``j < n`` (reedsolo convention)."""
    i = np.arange(p)[:, None]
    j = np.arange(n)[None, :]
    return GF_EXP[(i * (n - 1 - j)) % 255]


_GEN_CACHE: dict[tuple[int, int], np.ndarray] = {}


def _generator(k: int, p: int) -> np.ndarray:
    """Parity generator ``G`` (p x k) with ``parity = G @ data``."""
    key = (k, p)
    if key not in _GEN_CACHE:
        h = parity_check_matrix(k + p, p)
        _GEN_CACHE[key] = gf_matmul(gf_inv_matrix(h[:, k:]), h[:, :k])
    return _GEN_CACHE[key]


def rs_encode_matrix(data: np.ndarray, p: int) -> np.ndarray:
    """Return the ``(p, cols)`` parity rows for a ``(k, cols)`` data matrix."""
    if p == 0:
        return np.zeros((0, data.shape[1]), dtype=np.uint8)
    return gf_matmul(_generator(data.shape[0], p), data)


@dataclass
class BlockDecode:
    """Outcome of decoding one RS block."""

    ok: bool
    data: np.ndarray | None
    n_erasures: int
    corrected_rows: list[int] = field(default_factory=list)
    reason: str = ""


def _solve_erasures(h: np.ndarray, rows: np.ndarray, known: list[int], erased: list[int]) -> tuple[np.ndarray, bool]:
    """Solve for erased rows. Returns (filled matrix, consistent?)."""
    p = h.shape[0]
    e = len(erased)
    synd = gf_matmul(h[:, known], rows[known]) if known else np.zeros((p, rows.shape[1]), np.uint8)
    full = rows.copy()
    if e:
        x = gf_matmul(gf_inv_matrix(h[:e, erased]), synd[:e])
        full[erased] = x
        resid = gf_matmul(h[e:, erased], x) ^ synd[e:]
    else:
        resid = synd
    return full, not resid.any()


def rs_decode_matrix(rows: np.ndarray, present: np.ndarray, k: int, p: int) -> BlockDecode:
    """Decode one block.

    Parameters
    ----------
    rows: ``(k+p, cols)`` received symbols (erased rows may hold anything).
    present: boolean mask of rows that were received and passed the inner CRC.
    """
    n = k + p
    erased = [i for i in range(n) if not present[i]]
    known = [i for i in range(n) if present[i]]
    e = len(erased)
    if e > p:
        return BlockDecode(False, None, e, reason=f"{e} erasures > {p} parity")
    if p == 0:
        return BlockDecode(True, rows[:k].copy(), 0)
    h = parity_check_matrix(n, p)
    full, consistent = _solve_erasures(h, rows, known, erased)
    if consistent:
        return BlockDecode(True, full[:k], e)
    # Undetected bad row(s). Try exactly one extra bad row; need a spare syndrome to verify.
    if e + 1 >= p:
        return BlockDecode(False, None, e, reason="inconsistent syndromes, no spare parity to locate error")
    hits = []
    for r in known:
        kn = [i for i in known if i != r]
        cand, ok = _solve_erasures(h, rows, kn, sorted(erased + [r]))
        if ok:
            hits.append((r, cand))
    if len(hits) == 1:
        r, cand = hits[0]
        return BlockDecode(True, cand[:k], e, corrected_rows=[r])
    return BlockDecode(False, None, e, reason=f"inconsistent syndromes ({len(hits)} single-row candidates)")


# --------------------------------------------------------------------------- layout


def default_k_max(parity_permille: int) -> int:
    """Largest data-row count ``k`` with ``k + ceil(k * permille / 1000) <= 255``."""
    k = 255
    while k > 1 and k + parity_rows(k, parity_permille) > 255:
        k -= 1
    return k


def parity_rows(k: int, parity_permille: int) -> int:
    return (k * parity_permille + 999) // 1000


@dataclass(frozen=True)
class RSLayout:
    """Bit-packed arrangement of RS rows into oligos (DESIGN.md Part 0c).

    The stored bytes are cut into *rows* of ``row_bytes = ceil(payload_bits / 8)`` bytes and
    grouped into balanced RS blocks (``k`` data rows + ``p`` parity rows; each column of a
    block is one RS codeword). All rows (block0 data, block0 parity, block1 data, ...) are
    concatenated row-major into one bitstream, and oligo ``q`` carries bits
    ``[q * payload_bits, (q + 1) * payload_bits)``. Every oligo is therefore completely full,
    and density varies smoothly with the codec's capacity (no whole-byte quantisation).
    Because ``8 * row_bytes >= payload_bits``, a lost oligo erases at most one row in most
    columns of a block, plus at most two extra bytes where it straddles a row boundary.
    """

    stored_len: int
    payload_bits: int
    parity_permille: int
    k_max: int
    blocks: tuple[tuple[int, int], ...]

    @property
    def row_bytes(self) -> int:
        return -(-self.payload_bits // 8)

    @property
    def n_data_rows(self) -> int:
        return sum(k for k, _ in self.blocks)

    @property
    def n_parity_rows(self) -> int:
        return sum(p for _, p in self.blocks)

    @property
    def n_rows(self) -> int:
        return self.n_data_rows + self.n_parity_rows

    @property
    def stream_bits(self) -> int:
        return self.n_rows * self.row_bytes * 8

    @property
    def n_oligos(self) -> int:
        return -(-self.stream_bits // self.payload_bits)

    @property
    def data_oligo_equiv(self) -> float:
        """Oligos' worth of stream occupied by data rows (for net density)."""
        return self.n_data_rows * self.row_bytes * 8 / self.payload_bits

    def block_starts(self) -> list[int]:
        starts, q = [], 0
        for k, p in self.blocks:
            starts.append(q)
            q += k + p
        return starts


def plan_layout(stored_len: int, payload_bits: int, parity_permille: int, k_max: int | None = None) -> RSLayout:
    """Split ``stored_len`` bytes into RS rows/blocks for oligos carrying ``payload_bits`` each."""
    if not isinstance(stored_len, int) or stored_len < 0:
        raise ValueError("stored_len must be a nonnegative integer")
    if not isinstance(parity_permille, int) or not 0 <= parity_permille <= 65535:
        raise ValueError("parity_permille must be an integer in [0, 65535]")
    if payload_bits <= 0:
        raise ValueError("payload_bits must be positive")
    if k_max is None:
        k_max = default_k_max(parity_permille)
    if not isinstance(k_max, int) or not 1 <= k_max <= 255:
        raise ValueError("k_max must be an integer in [1, 255]")
    if k_max + parity_rows(k_max, parity_permille) > 255:
        raise ValueError("k_max too large for GF(256) RS with this parity ratio")
    rb = -(-payload_bits // 8)
    n_data = -(-stored_len // rb)
    blocks: list[tuple[int, int]] = []
    if n_data:
        nb = -(-n_data // k_max)
        base, extra = divmod(n_data, nb)
        for i in range(nb):
            k = base + (1 if i < extra else 0)
            blocks.append((k, parity_rows(k, parity_permille)))
    return RSLayout(stored_len, payload_bits, parity_permille, k_max, tuple(blocks))


def _bits_to_int(bits: np.ndarray) -> int:
    n = len(bits)
    pad = (-n) % 8
    return int.from_bytes(np.packbits(bits).tobytes(), "big") >> pad if n else 0


def _int_to_bits(v: int, n: int) -> np.ndarray:
    nb = -(-n // 8)
    pad = nb * 8 - n
    return np.unpackbits(np.frombuffer((v << pad).to_bytes(nb, "big"), dtype=np.uint8))[:n]


def rs_encode(stored: bytes, layout: RSLayout) -> list[int]:
    """Return each oligo's ``payload_bits``-bit chunk of the RS-coded bitstream, as ints."""
    rb = layout.row_bytes
    padded = stored + bytes(layout.n_data_rows * rb - len(stored))
    mat = np.frombuffer(padded, dtype=np.uint8).reshape(layout.n_data_rows, rb) if layout.n_data_rows else np.zeros((0, rb), np.uint8)
    rows: list[np.ndarray] = []
    r = 0
    for k, p in layout.blocks:
        data = mat[r : r + k]
        r += k
        rows.append(data)
        rows.append(rs_encode_matrix(data, p))
    if not rows:
        return []
    bits = np.unpackbits(np.concatenate(rows).reshape(-1))
    d = layout.payload_bits
    bits = np.concatenate([bits, np.zeros(layout.n_oligos * d - len(bits), np.uint8)])
    return [_bits_to_int(bits[q * d : (q + 1) * d]) for q in range(layout.n_oligos)]


def _ee_column(word: np.ndarray, erased: list[int], p: int) -> np.ndarray | None:
    """Errors-and-erasures decode of one codeword (reedsolo, same convention as ours)."""
    import reedsolo

    if reedsolo.gf_exp[1] != 2 or len(reedsolo.gf_exp) < 255:
        reedsolo.init_tables(prim=_PRIM, generator=2, c_exp=8)
    try:
        msg, _, _ = reedsolo.rs_correct_msg(bytearray(word.tobytes()), p, fcr=0, generator=2,
                                            erase_pos=erased or None)
    except reedsolo.ReedSolomonError:
        return None
    full = np.frombuffer(bytes(msg), dtype=np.uint8)
    return np.concatenate([full, rs_encode_matrix(full[:, None], p)[:, 0]]) if p else full


def _decode_block_masked(h: np.ndarray, rows: np.ndarray, known: np.ndarray, k: int, p: int,
                         errors: bool = False) -> tuple[np.ndarray | None, str, int]:
    """Erasure-decode one block with a per-byte known mask. Columns are grouped by erasure
    pattern so each distinct pattern needs one solve. Returns (data rows, "", n_ee_columns)
    or (None, reason, n). An inconsistent spare syndrome is reported as ``"inconsistent"``,
    unless ``errors=True``: then each column of an inconsistent pattern is decoded with
    errors and erasures (2 * errors + erasures <= p), and a column that fails is reported."""
    n = k + p
    out = rows.copy()
    n_ee = 0
    patterns: dict[tuple[int, ...], list[int]] = {}
    for c in range(rows.shape[1]):
        patterns.setdefault(tuple(np.nonzero(~known[:, c])[0].tolist()), []).append(c)
    for erased, cols in patterns.items():
        if len(erased) > p:
            return None, f"{len(erased)} erasures > {p} parity", n_ee
        if p == 0:
            continue
        kn = [i for i in range(n) if i not in erased]
        sub = rows[:, cols]
        full, ok = _solve_erasures(h, sub, kn, list(erased))
        if ok:
            out[:, cols] = full
            continue
        if not errors:
            return None, "inconsistent", n_ee
        for c in cols:
            # cheap per-column consistency check first (pattern group may mix good columns)
            fc, okc = _solve_erasures(h, rows[:, [c]], kn, list(erased))
            if okc:
                out[:, c] = fc[:, 0]
                continue
            n_ee += 1
            word = _ee_column(rows[:, c], list(erased), p)
            if word is None:
                return None, f"errors-and-erasures failed ({len(erased)} erasures)", n_ee
            out[:, c] = word
    return out[:k], "", n_ee


@dataclass
class RSDecodeResult:
    data: bytes | None
    blocks_failed: int
    erasures_per_block: list[int]  # max erased rows over the block's columns
    corrected_rows: int  # undetected-bad oligos located and removed
    reasons: list[str]
    ee_columns: int = 0  # columns that needed errors-and-erasures decoding


def rs_decode(received: dict[int, int], layout: RSLayout, unverified: dict[int, int] | None = None,
              errors: bool = False) -> RSDecodeResult:
    """Decode from ``{oligo position: payload_bits int}`` (CRC-valid oligos).

    ``unverified`` holds payload bits of oligos that failed the inner CRC (best decoded
    guess). With ``errors=True`` they fill their rows as *possibly wrong* symbols instead of
    erasures, and inconsistent columns are decoded with errors and erasures. A substitution
    that corrupts a few bytes of an oligo then costs 2 parity symbols in only the affected
    columns, rather than 1 erasure in every column. Without ``errors``, ``unverified`` is
    ignored (erasure-only decoding, plus the single-undetected-bad-oligo search).
    """
    d, rb = layout.payload_bits, layout.row_bytes
    total = layout.n_oligos * d
    bits = np.zeros(total, np.uint8)
    have = np.zeros(total, bool)
    filled = dict(received)
    if errors and unverified:
        for q, v in unverified.items():
            filled.setdefault(q, v)
    for q, v in filled.items():
        if 0 <= q < layout.n_oligos:
            bits[q * d : (q + 1) * d] = _int_to_bits(v, d)
            have[q * d : (q + 1) * d] = True
    sb = layout.stream_bits
    stream = np.packbits(bits[:sb]).reshape(layout.n_rows, rb) if layout.n_rows else np.zeros((0, rb), np.uint8)
    known = have[:sb].reshape(-1, 8).all(axis=1).reshape(layout.n_rows, rb) if layout.n_rows else np.zeros((0, rb), bool)
    known_v = known
    if errors and unverified:  # mask of CRC-verified bytes only, for the erasure fallback
        hv = np.zeros(total, bool)
        for q in received:
            if 0 <= q < layout.n_oligos:
                hv[q * d : (q + 1) * d] = True
        known_v = hv[:sb].reshape(-1, 8).all(axis=1).reshape(layout.n_rows, rb)
    chunks, failed, erasures, corrected, reasons, n_ee = [], 0, [], 0, [], 0
    for start, (k, p) in zip(layout.block_starts(), layout.blocks):
        n = k + p
        rows, kn = stream[start : start + n], known[start : start + n]
        erasures.append(int((~kn).sum(axis=0).max()) if n else 0)
        h = parity_check_matrix(n, p) if p else None
        data, reason, ne = _decode_block_masked(h, rows, kn, k, p, errors)
        n_ee += ne
        if data is None and known_v is not known:
            # too many wrong unverified rows: they cost 2 parity each as errors, 1 as erasures
            data, reason, _ = _decode_block_masked(h, rows, known_v[start : start + n], k, p)
        if data is None and reason == "inconsistent":
            # one undetected bad oligo: try erasing each received oligo overlapping this block
            b0, b1 = start * rb * 8, (start + n) * rb * 8
            hits = []
            for q in sorted(received):
                lo, hi = q * d, (q + 1) * d
                if hi <= b0 or lo >= b1:
                    continue
                h2 = have.copy()
                h2[lo:hi] = False
                kn2 = h2[:sb].reshape(-1, 8).all(axis=1).reshape(layout.n_rows, rb)[start : start + n]
                cand, why, _ = _decode_block_masked(h, rows, kn2, k, p)
                if cand is not None:
                    hits.append(cand)
            if len(hits) == 1:
                data, corrected = hits[0], corrected + 1
            else:
                reason = f"inconsistent syndromes ({len(hits)} single-oligo candidates)"
        if data is None:
            failed += 1
            reasons.append(reason)
            chunks.append(bytes(k * rb))
        else:
            chunks.append(data.tobytes())
    if failed:
        return RSDecodeResult(None, failed, erasures, corrected, reasons, n_ee)
    return RSDecodeResult(b"".join(chunks)[: layout.stored_len], 0, erasures, corrected, reasons, n_ee)


# --------------------------------------------------------------------------- inner CRC

CRC_BYTES = 4
CRC_BITS = 32
_CRC_MASK = (1 << CRC_BITS) - 1


def crc(data: bytes) -> int:
    """CRC-32 (IEEE, ``zlib.crc32``)."""
    return zlib.crc32(data) & _CRC_MASK


def add_crc(index: int, payload: bytes) -> bytes:
    """``payload || crc(index_be32 || payload)``."""
    return payload + crc(index.to_bytes(4, "big") + payload).to_bytes(CRC_BYTES, "big")


def check_crc(index: int, frame: bytes) -> bytes | None:
    """Return the payload if the frame's CRC matches, else ``None``."""
    if len(frame) < CRC_BYTES:
        return None
    payload, tag = frame[:-CRC_BYTES], frame[-CRC_BYTES:]
    if crc(index.to_bytes(4, "big") + payload).to_bytes(CRC_BYTES, "big") != tag:
        return None
    return payload


def add_crc_bits(index: int, data: int, nbits: int) -> int:
    """Frame = ``data << CRC_BITS | crc(index_be32 || data as ceil(nbits/8) big-endian bytes)``."""
    tag = crc(index.to_bytes(4, "big") + data.to_bytes(-(-nbits // 8), "big"))
    return (data << CRC_BITS) | tag


def check_crc_bits(index: int, frame: int, frame_bits: int) -> int | None:
    """Return the ``frame_bits - CRC_BITS`` data bits if the CRC matches, else ``None``."""
    nbits = frame_bits - CRC_BITS
    if frame < 0 or frame >> frame_bits:
        return None
    data, tag = frame >> CRC_BITS, frame & _CRC_MASK
    if crc(index.to_bytes(4, "big") + data.to_bytes(-(-nbits // 8), "big")) != tag:
        return None
    return data
