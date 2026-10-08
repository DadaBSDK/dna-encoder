from __future__ import annotations

import os

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dnastore.ecc import (
    add_crc,
    check_crc,
    default_k_max,
    gf_inv_matrix,
    gf_matmul,
    parity_check_matrix,
    plan_layout,
    rs_decode,
    rs_decode_matrix,
    rs_encode,
    rs_encode_matrix,
)


def test_gf_inverse():
    rng = np.random.default_rng(0)
    for n in (1, 3, 8):
        m = parity_check_matrix(n + 5, n)[:, 2 : 2 + n]  # Vandermonde, invertible
        assert np.array_equal(gf_matmul(m, gf_inv_matrix(m)), np.eye(n, dtype=np.uint8))
    del rng


@given(st.integers(1, 40), st.integers(0, 12), st.data())
@settings(max_examples=60, deadline=None)
def test_rs_any_p_erasures_recoverable(k, p, data):
    rng = np.random.default_rng(data.draw(st.integers(0, 2**32 - 1)))
    msg = rng.integers(0, 256, size=(k, 5), dtype=np.uint8)
    code = np.concatenate([msg, rs_encode_matrix(msg, p)])
    assert not gf_matmul(parity_check_matrix(k + p, p), code).any() if p else True
    e = data.draw(st.integers(0, p))
    erased = rng.choice(k + p, size=e, replace=False)
    present = np.ones(k + p, bool)
    present[erased] = False
    rx = code.copy()
    rx[erased] = rng.integers(0, 256, size=(e, 5), dtype=np.uint8)
    res = rs_decode_matrix(rx, present, k, p)
    assert res.ok and np.array_equal(res.data, msg)


def test_rs_too_many_erasures_fails():
    msg = np.arange(20, dtype=np.uint8).reshape(10, 2)
    code = np.concatenate([msg, rs_encode_matrix(msg, 3)])
    present = np.ones(13, bool)
    present[:4] = False
    assert not rs_decode_matrix(code, present, 10, 3).ok


def test_rs_single_undetected_row_error_corrected():
    rng = np.random.default_rng(1)
    msg = rng.integers(0, 256, size=(30, 8), dtype=np.uint8)
    code = np.concatenate([msg, rs_encode_matrix(msg, 6)])
    code[7, 3] ^= 0x5A  # corrupted row that "passed" CRC
    present = np.ones(36, bool)
    present[[1, 2]] = False
    res = rs_decode_matrix(code, present, 30, 6)
    assert res.ok and res.corrected_rows == [7] and np.array_equal(res.data, msg)


def test_layout_balanced_and_bounded():
    lay = plan_layout(100_000, 203, 150)
    assert all(k + p <= 255 for k, p in lay.blocks)
    ks = [k for k, _ in lay.blocks]
    assert max(ks) - min(ks) <= 1
    rb = lay.row_bytes
    assert rb == 26 and lay.n_data_rows * rb >= 100_000 > (lay.n_data_rows - 1) * rb
    assert lay.n_oligos == -(-lay.stream_bits // 203)
    assert default_k_max(150) + -(-default_k_max(150) * 150 // 1000) <= 255
    assert plan_layout(0, 203, 150).n_oligos == 0


@pytest.mark.parametrize("size,bits", [(1, 203), (25, 200), (26, 203), (5000, 203), (5000, 97), (7777, 251)])
def test_rs_bitpacked_roundtrip_with_dropout(size, bits):
    """Bit-packed stream: any oligo payload width; dropping oligos within parity recovers."""
    stored = os.urandom(size)
    lay = plan_layout(size, bits, 200)
    chunks = rs_encode(stored, lay)
    assert all(0 <= c < 2**bits for c in chunks) and len(chunks) == lay.n_oligos
    rx = dict(enumerate(chunks))
    assert rs_decode(rx, lay).data == stored
    rng = np.random.default_rng(size)
    # drop a fraction well inside the parity budget (each lost oligo can straddle 2 rows)
    p_min = min(p for _, p in lay.blocks)
    for q in rng.choice(lay.n_oligos, size=max(0, min(lay.n_oligos - 1, p_min // 2 * len(lay.blocks))), replace=False):
        rx.pop(int(q))
    res = rs_decode(rx, lay)
    assert res.data == stored, res.reasons


def test_rs_bitpacked_locates_one_undetected_bad_oligo():
    stored = os.urandom(3000)
    lay = plan_layout(len(stored), 203, 200)
    chunks = rs_encode(stored, lay)
    rx = dict(enumerate(chunks))
    rx[5] ^= (1 << 100) | 1  # corrupted oligo that "passed" CRC
    res = rs_decode(rx, lay)
    assert res.data == stored and res.corrected_rows == 1


def test_bitpacking_smooth_in_payload_bits():
    """Oligo count decreases monotonically (no whole-byte plateaus) as capacity grows."""
    counts = [plan_layout(20_000, b, 150).n_oligos for b in range(180, 260, 4)]
    assert all(a >= b for a, b in zip(counts, counts[1:]))
    assert len(set(counts)) > 12


def test_crc_binds_index():
    frame = add_crc(1234, b"payload")
    assert check_crc(1234, frame) == b"payload"
    assert check_crc(1235, frame) is None
    bad = bytearray(frame)
    bad[0] ^= 1
    assert check_crc(1234, bytes(bad)) is None


@given(st.integers(1, 200), st.integers(1, 54), st.data())
@settings(max_examples=80, deadline=None)
def test_rs_codewords_identical_to_reedsolo(k, p, data):
    """The custom RS must produce exactly reedsolo's systematic codewords (fcr=0, 0x11d, alpha=2)."""
    reedsolo = pytest.importorskip("reedsolo")
    if k + p > 255:
        k = 255 - p
    rng = np.random.default_rng(data.draw(st.integers(0, 2**32 - 1)))
    msg = rng.integers(0, 256, size=(k, 3), dtype=np.uint8)
    ours = rs_encode_matrix(msg, p)
    rs = reedsolo.RSCodec(p, nsize=255, fcr=0, prim=0x11D, generator=2, c_exp=8)
    for col in range(3):
        ref = bytes(rs.encode(bytes(msg[:, col])))
        assert ref == bytes(msg[:, col]) + bytes(ours[:, col])
