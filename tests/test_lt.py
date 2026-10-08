"""Tests for the LT fountain core (dnastore/lt.py)."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dnastore.lt import (
    LTDecoder,
    RobustSoliton,
    droplet_neighbors,
    encode_droplet,
    generate_screened,
)


def _segments(k: int, seg_len: int, seed: int) -> list[bytes]:
    rng = np.random.default_rng(seed)
    return [rng.bytes(seg_len) for _ in range(k)]


def test_robust_soliton_is_distribution():
    for k in (1, 2, 10, 1000):
        d = RobustSoliton(k)
        assert d.pmf.shape == (k,)
        assert np.all(d.pmf >= 0)
        assert abs(d.pmf.sum() - 1.0) < 1e-12
        assert d.cdf[-1] == 1.0


def test_neighbors_deterministic_and_valid():
    d = RobustSoliton(500)
    for seed in (0, 1, 2**32 - 1, 123456789):
        n1 = droplet_neighbors(seed, 500, d)
        assert n1 == droplet_neighbors(seed, 500, d)
        assert len(set(n1)) == len(n1) >= 1
        assert all(0 <= j < 500 for j in n1) and n1 == sorted(n1)
    with pytest.raises(ValueError):
        droplet_neighbors(2**32, 500, d)


def test_degree_histogram_matches_pmf():
    k = 200
    d = RobustSoliton(k)
    degs = np.array([len(droplet_neighbors(s, k, d)) for s in range(20000)])
    assert abs(degs.mean() - d.mean_degree) < 0.05 * d.mean_degree
    assert abs((degs == 1).mean() - d.pmf[0]) < 0.01


def test_encode_droplet_is_xor_of_neighbors():
    segs = _segments(50, 16, 1)
    d = RobustSoliton(50)
    nb = droplet_neighbors(77, 50, d)
    ref = np.zeros(16, np.uint8)
    for j in nb:
        ref ^= np.frombuffer(segs[j], np.uint8)
    assert encode_droplet(segs, 77, d) == ref.tobytes()


@pytest.mark.parametrize("k,overhead", [(50, 0.25), (300, 0.12), (2000, 0.08)])
def test_decode_with_overhead(k, overhead):
    segs = _segments(k, 20, k)
    d = RobustSoliton(k)
    dec = LTDecoder(k, 20, d)
    for s in range(int(k * (1 + overhead))):
        dec.add(s, encode_droplet(segs, s, d))
    r = dec.decode()
    assert r.ok and r.n_missing == 0 and r.segments == segs


@given(st.integers(20, 120), st.integers(0, 2**31))
@settings(max_examples=25, deadline=None)
def test_decode_random_property(k, base):
    """With generous overhead decoding succeeds; any output is always correct."""
    segs = _segments(k, 8, base)
    d = RobustSoliton(k)
    dec = LTDecoder(k, 8, d)
    for s in range(base, base + int(k * 1.5)):
        dec.add(s % 2**32, encode_droplet(segs, s % 2**32, d))
    r = dec.decode()
    assert r.ok and r.segments == segs


def test_too_few_droplets_reported_not_silent():
    k = 400
    segs = _segments(k, 10, 3)
    d = RobustSoliton(k)
    dec = LTDecoder(k, 10, d)
    for s in range(k // 2):
        dec.add(s, encode_droplet(segs, s, d))
    r = dec.decode()
    assert not r.ok and r.segments is None and r.n_missing > 0


def test_duplicates_ignored_and_counted():
    k = 40
    segs = _segments(k, 6, 4)
    d = RobustSoliton(k)
    dec = LTDecoder(k, 6, d)
    for s in range(80):
        x = encode_droplet(segs, s, d)
        dec.add(s, x)
        dec.add(s, x)
    r = dec.decode()
    assert r.n_duplicates == 80 and r.n_droplets == 80
    assert r.ok and r.segments == segs


def test_gaussian_elimination_fallback_when_peeling_stalls():
    """No degree-1 droplets at all: peeling recovers nothing, GE must solve everything."""
    k = 8
    segs = _segments(k, 12, 5)
    d = RobustSoliton(k)
    chosen, rows = [], []
    rank_mat = np.zeros((0, k), dtype=np.uint8)
    s = 0
    while len(chosen) < 3 * k:
        nb = droplet_neighbors(s, k, d)
        if len(nb) >= 2:
            chosen.append(s)
            v = np.zeros(k, np.uint8)
            v[nb] = 1
            rows.append(v)
        s += 1
    # confirm full GF(2) rank of the chosen system
    m = np.array(rows) % 2
    r = 0
    for c in range(k):
        piv = next((i for i in range(r, len(m)) if m[i, c]), None)
        if piv is None:
            continue
        m[[r, piv]] = m[[piv, r]]
        for i in range(len(m)):
            if i != r and m[i, c]:
                m[i] ^= m[r]
        r += 1
    assert r == k
    dec = LTDecoder(k, 12, d)
    for s in chosen:
        dec.add(s, encode_droplet(segs, s, d))
    res = dec.decode()
    assert res.n_peeled == 0 and res.n_ge == k
    assert res.ok and res.segments == segs


def test_partial_ge_reports_remaining():
    """Rank-deficient residual: GE solves what is determined and reports the rest missing."""
    k = 6
    segs = _segments(k, 4, 6)
    d = RobustSoliton(k)
    # droplets only over segments {0..3}, so segments 4 and 5 can never be solved
    chosen = [s for s in range(5000) if set(droplet_neighbors(s, k, d)) <= {0, 1, 2, 3}][:20]
    dec = LTDecoder(k, 4, d)
    for s in chosen:
        dec.add(s, encode_droplet(segs, s, d))
    r = dec.decode()
    assert not r.ok and r.n_missing >= 2


def test_screening_hook_counts_rejections():
    k = 30
    segs = _segments(k, 10, 7)
    d = RobustSoliton(k)
    stats: dict = {}
    accept = lambda seed, data: data[0] % 2 == 0  # noqa: E731  (rejects about half)
    out = generate_screened(segs, iter(range(10_000)), accept, 60, d, stats)
    assert len(out) == 60
    assert all(x[1][0] % 2 == 0 for x in out)
    assert stats["accepted"] == 60 and stats["rejected"] > 0
    assert stats["drawn"] == stats["accepted"] + stats["rejected"]
    for seed, data in out:
        assert encode_droplet(segs, seed, d) == data
    # droplets screened by a rule that is not linear in the data still decode
    nonlinear = lambda seed, data: (data[0] * 7 + data[1]) % 5 != 0  # noqa: E731
    dec = LTDecoder(k, 10, d)
    for seed, data in generate_screened(segs, iter(range(10_000, 20_000)), nonlinear, 3 * k, d):
        dec.add(seed, data)
    assert dec.decode().segments == segs


def test_linear_screening_rule_makes_system_rank_deficient():
    """Pitfall (documented in docs/notes/fountain.md): an accept rule that is a GF(2)-linear
    function of the payload (here: parity of one bit) only admits droplets whose neighbour
    sets have even overlap with the set of segments having that bit set. All accepted rows
    are then orthogonal to that indicator vector, so the segments cannot all be recovered,
    however many droplets are collected."""
    k = 30
    segs = _segments(k, 10, 7)
    d = RobustSoliton(k)
    parity = lambda seed, data: data[0] % 2 == 0  # noqa: E731
    dec = LTDecoder(k, 10, d)
    for seed, data in generate_screened(segs, iter(range(50_000)), parity, 10 * k, d):
        dec.add(seed, data)
    r = dec.decode()
    odd = sum(sg[0] & 1 for sg in segs)
    if odd:  # if no segment has the bit set, the constraint is vacuous
        assert not r.ok


def test_screening_hook_raises_when_seeds_run_out():
    segs = _segments(10, 4, 8)
    with pytest.raises(RuntimeError):
        generate_screened(segs, iter(range(5)), lambda s, x: True, 10)
