from __future__ import annotations

import itertools

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dnastore.dnautil import revcomp
from dnastore.kmer import (
    KmerModel,
    adjacent_contrast,
    code_to_kmer,
    confusability_table,
    encode_kmers,
    kmer_code,
    load_model,
    oligo_confusability,
    rc_code_table,
)

dna = st.text(alphabet="ACGT", min_size=0, max_size=60)


def _fake(k: int, levels, std=None, units="std") -> KmerModel:
    return KmerModel("fake", k, np.asarray(levels, float), None if std is None else np.asarray(std, float), units, "x")


def _brute_conf(model, thr_fn):
    k = model.k
    out = []
    for c in range(4**k):
        x = code_to_kmer(c, k)
        n = 0
        for pos in range(k):
            for b in "ACGT":
                if b != x[pos]:
                    y = x[:pos] + b + x[pos + 1 :]
                    d = kmer_code(y)
                    n += abs(model.levels[c] - model.levels[d]) < thr_fn(c, d)
        out.append(n / (3 * k))
    return np.array(out)


def test_conf_hand_example_k1():
    # k=1: A=0.0, C=0.04, G=1.0, T=1.03; delta=0.05 -> A~C, G~T only.
    m = _fake(1, [0.0, 0.04, 1.0, 1.03])
    t = confusability_table(m, delta=0.05, use_cache=False)
    assert np.allclose(t, [1 / 3, 1 / 3, 1 / 3, 1 / 3])
    t = confusability_table(m, delta=2.0, use_cache=False)
    assert np.allclose(t, 1.0)


@pytest.mark.parametrize("k", [2, 3])
def test_conf_matches_bruteforce_std(k):
    rng = np.random.default_rng(k)
    m = _fake(k, rng.normal(size=4**k))
    t = confusability_table(m, delta=0.3, use_cache=False)
    assert np.allclose(t, _brute_conf(m, lambda c, d: 0.3))
    assert t.shape == (4**k,) and t.min() >= 0 and t.max() <= 1


def test_conf_matches_bruteforce_pa_pooled_std():
    rng = np.random.default_rng(9)
    lv, sd = rng.normal(80, 10, 64), rng.uniform(1, 3, 64)
    m = _fake(3, lv, sd, units="pA")
    t = confusability_table(m, c=1.5, use_cache=False)
    ref = _brute_conf(m, lambda a, b: 1.5 * np.sqrt((sd[a] ** 2 + sd[b] ** 2) / 2))
    assert np.allclose(t, ref)


def test_conf_param_units_enforced():
    with pytest.raises(ValueError):
        confusability_table(_fake(1, [0, 1, 2, 3]), c=1.0, use_cache=False)
    with pytest.raises(ValueError):
        confusability_table(_fake(1, [0, 1, 2, 3], [1, 1, 1, 1], "pA"), delta=1.0, use_cache=False)


@given(dna, st.integers(1, 6))
def test_rolling_codes_match_direct(seq, k):
    codes = encode_kmers(seq, k)
    assert list(codes) == [kmer_code(seq[i : i + k]) for i in range(len(seq) - k + 1)]


@pytest.mark.parametrize("k", [1, 2, 3, 5])
def test_rc_code_table(k):
    rc = rc_code_table(k)
    for c in range(min(4**k, 300)):
        assert code_to_kmer(int(rc[c]), k) == revcomp(code_to_kmer(c, k))
    assert np.array_equal(rc[rc], np.arange(4**k))


@given(st.text(alphabet="ACGT", min_size=12, max_size=60))
@settings(max_examples=50)
def test_both_strands_via_rc_table(seq):
    k = 3
    rng = np.random.default_rng(0)
    table = rng.random(4**k).astype(np.float32)
    r = oligo_confusability(seq, table, k)
    rc_direct = table[encode_kmers(revcomp(seq), k)]
    assert np.isclose(r["rc_mean"], rc_direct.mean()) and np.isclose(r["rc_max"], rc_direct.max())
    # symmetric under reverse complement of the input
    r2 = oligo_confusability(revcomp(seq), table, k)
    assert np.isclose(r["mean"], r2["mean"]) and np.isclose(r["max"], r2["max"])


def test_adjacent_contrast_hand_example():
    m = _fake(1, [0.0, 1.0, 3.0, 6.0], units="pA")  # A C G T
    r = adjacent_contrast("ACGT", m, tau=1.5)
    # fwd diffs: |0-1|=1, |1-3|=2, |3-6|=3 ; rc = ACGT as well
    assert r["fwd_min"] == 1.0 and np.isclose(r["fwd_mean"], 2.0)
    assert np.isclose(r["fwd_frac_below_tau"], 1 / 3)
    r = adjacent_contrast("AAC", m, tau=1.0)  # diffs 0, 1 -> both <= tau
    assert r["fwd_frac_below_tau"] == 1.0
    rc = adjacent_contrast(revcomp("AAC"), m, tau=1.0)
    assert np.isclose(r["rc_mean"], rc["fwd_mean"])


def test_real_models_load_and_cover():
    for name, k, units in (("r10.4.1_9mer", 9, "std"), ("r9.4.1_6mer", 6, "pA")):
        try:
            m = load_model(name)
        except FileNotFoundError:
            pytest.skip("k-mer tables not downloaded")
        assert m.k == k and m.units == units and m.levels.shape == (4**k,)
        assert np.isfinite(m.levels).all()
    r9 = load_model("r9.4.1_6mer")
    assert r9.levels[kmer_code("AAAAAA")] == pytest.approx(86.486336)
    assert r9.level_std[kmer_code("AAAAAC")] == pytest.approx(1.517846)
    r10 = load_model("r10.4.1_9mer")
    assert r10.levels[kmer_code("AAAAAAAAC")] == pytest.approx(-1.6519798040390015)


def test_real_table_spotcheck_r10():
    try:
        m = load_model("r10.4.1_9mer")
    except FileNotFoundError:
        pytest.skip("k-mer tables not downloaded")
    t = confusability_table(m, delta=0.05)
    rng = np.random.default_rng(1)
    for c in rng.integers(0, 4**9, 20):
        x = code_to_kmer(int(c), 9)
        n = sum(
            abs(m.levels[c] - m.levels[kmer_code(x[:p] + b + x[p + 1 :])]) < 0.05
            for p, b in itertools.product(range(9), "ACGT")
            if b != x[p]
        )
        assert np.isclose(t[c], n / 27)
