from __future__ import annotations

import numpy as np
import pytest

from dnastore.channel import ChannelParams, conf_multipliers, simulate

RNG = np.random.default_rng(123)


def _oligos(n=200, L=200, seed=0):
    r = np.random.default_rng(seed)
    return ["".join("ACGT"[i] for i in r.integers(0, 4, L)) for _ in range(n)]


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_error_rates_match_params():
    ol = _oligos()
    p = ChannelParams(p_sub=0.02, p_ins=0.01, p_del=0.015, mean_coverage=20, seed=1)
    reads, truth = simulate(ol, p)
    nb = sum(len(ol[r.oligo]) for r in reads)
    for attr, rate in (("n_sub", 0.02), ("n_ins", 0.01), ("n_del", 0.015)):
        k = sum(getattr(r, attr) for r in reads)
        sd = np.sqrt(nb * rate * (1 - rate))
        assert abs(k - nb * rate) < 5 * sd, attr
    # read lengths consistent with recorded events
    for r in reads[:500]:
        assert len(r.seq) == len(ol[r.oligo]) + r.n_ins - r.n_del


def test_noiseless_reads_equal_oligo_or_rc_and_rc_fraction():
    ol = _oligos(100)
    reads, _ = simulate(ol, ChannelParams(mean_coverage=40, seed=2))
    assert len(reads) == 4000
    for r in reads:
        assert r.seq == (_rc(ol[r.oligo]) if r.reverse else ol[r.oligo])
    frac = np.mean([r.reverse for r in reads])
    assert abs(frac - 0.5) < 5 * np.sqrt(0.25 / len(reads))


def test_reproducible_and_seed_sensitive():
    ol = _oligos(50)
    p = ChannelParams(p_sub=0.01, p_ins=0.01, p_del=0.01, seed=7)
    a, _ = simulate(ol, p)
    b, _ = simulate(ol, p)
    c, _ = simulate(ol, ChannelParams(p_sub=0.01, p_ins=0.01, p_del=0.01, seed=8))
    assert [r.seq for r in a] == [r.seq for r in b]
    assert [r.seq for r in a] != [r.seq for r in c]


def test_dropout_and_depth():
    ol = _oligos(2000, L=50)
    reads, truth = simulate(ol, ChannelParams(dropout=0.2, mean_coverage=5, seed=3))
    assert abs(truth.dropped.mean() - 0.2) < 5 * np.sqrt(0.16 / 2000)
    assert all(not truth.dropped[r.oligo] for r in reads)
    assert len(reads) == 10000  # fixed total depth, redistributed to survivors


def test_lognormal_and_model_coverage():
    ol = _oligos(3000, L=60)
    _, t1 = simulate(ol, ChannelParams(coverage="lognormal", lognormal_sigma=0.8, mean_coverage=30, seed=4))
    _, t0 = simulate(ol, ChannelParams(coverage="poisson", mean_coverage=30, seed=4))
    assert np.std(t1.read_counts) > 1.5 * np.std(t0.read_counts)
    gc = np.array([(s.count("G") + s.count("C")) / len(s) for s in ol])
    acc = np.random.default_rng(0).uniform(0, 1, len(ol))
    p = ChannelParams(coverage="model", coverage_preset="acc_assumed_strong", mean_coverage=30, seed=5)
    _, tm = simulate(ol, p, features={"gc": gc, "accessibility": acc})
    _, tz = simulate(ol, ChannelParams(**{**p.__dict__, "zero_accessibility": True}), features={"gc": gc, "accessibility": acc})
    lo = acc < 0.2
    # the assumed knob depresses low-accessibility oligos; zeroing it removes the effect
    assert tm.weights[lo].mean() / tm.weights[~lo].mean() < 0.9
    assert 0.8 < tz.weights[lo].mean() / tz.weights[~lo].mean() < 1.25


def test_conf_weighted_mean_rate_and_targeting():
    pytest.importorskip("dnastore.kmer")
    ol = _oligos(300)
    p = ChannelParams(p_sub=0.02, mean_coverage=10, conf_weighted=True, seed=9)
    reads, _ = simulate(ol, p)
    nb = sum(len(ol[r.oligo]) for r in reads)
    rate = sum(r.n_sub for r in reads) / nb
    assert abs(rate - 0.02) < 0.002  # normalised: mean rate ~ base rate on random sequences
    m = np.concatenate([conf_multipliers(s, p) for s in ol])
    assert abs(m.mean() - 1.0) < 0.02
    assert m.max() > 1.5 and m.min() < 0.7  # it does redistribute errors toward confusable k-mers
