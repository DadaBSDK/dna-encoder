from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from dnastore.dnautil import (
    bases_to_bytes_2bit,
    bytes_to_bases_2bit,
    bytes_to_trits,
    keystream,
    max_bytes_for_trits,
    min_trits_for_bytes,
    revcomp,
    rotate_decode,
    rotate_encode,
    trits_to_bytes,
)


def test_revcomp():
    assert revcomp("AACGT") == "ACGTT"
    assert revcomp(revcomp("GATTACA")) == "GATTACA"


@given(st.binary(max_size=64))
def test_2bit_roundtrip(data):
    assert bases_to_bytes_2bit(bytes_to_bases_2bit(data)) == data


def test_2bit_mapping():
    assert bytes_to_bases_2bit(bytes([0b00011011])) == "ACGT"


@given(st.integers(1, 300))
def test_max_bytes_for_trits_is_tight(n):
    nb = max_bytes_for_trits(n)
    assert 256**nb <= 3**n < 256 ** (nb + 1)


@given(st.binary(max_size=40))
def test_trits_roundtrip(data):
    n = min_trits_for_bytes(len(data))
    assert trits_to_bytes(bytes_to_trits(data, n), len(data)) == data


def test_trit_overflow_detected():
    n = 6  # 3**6 = 729 > 256: values >= 256 are invalid for 1 byte
    assert trits_to_bytes([2] * n, 1) is None


@given(st.lists(st.integers(0, 2), max_size=200), st.sampled_from("ACGT"))
@settings(max_examples=50)
def test_rotate_roundtrip_and_no_repeats(trits, prev):
    seq = rotate_encode(trits, prev)
    full = prev + seq
    assert all(a != b for a, b in zip(full, full[1:]))
    assert rotate_decode(seq, prev) == trits


def test_rotate_decode_rejects_repeat():
    assert rotate_decode("AA", "C") is None


def test_keystream_deterministic_and_keyed():
    assert keystream(16, 1, 2) == keystream(16, 1, 2)
    assert keystream(16, 1, 2) != keystream(16, 1, 3)
    assert keystream(0, 1) == b""
