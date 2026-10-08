from __future__ import annotations

import edlib
import numpy as np
import pytest

from dnastore.synthesis import PRESETS, ZERO, SynthesisModel, apply_synthesis_errors, mutate


def _seqs(n, L, seed=0):
    rng = np.random.default_rng(seed)
    return ["".join("ACGT"[x] for x in rng.integers(0, 4, L)) for _ in range(n)]


def test_zero_model_is_identity():
    s = _seqs(20, 150)
    out = apply_synthesis_errors(s, ZERO, np.random.default_rng(1), copies=3)
    assert all(m == o for ms, o in zip(out, s) for m in ms)


def test_reproducible_with_seed():
    s = _seqs(10, 150)
    m = PRESETS["photolitho_gimpel2024"]
    a = apply_synthesis_errors(s, m, np.random.default_rng(7), copies=2)
    b = apply_synthesis_errors(s, m, np.random.default_rng(7), copies=2)
    assert a == b


def test_every_preset_has_source_and_valid_rates():
    for name, m in PRESETS.items():
        assert m.name == name and m.source
        assert 0 <= m.p_del < 1 and 0 <= m.p_ins < 1 and 0 <= m.p_sub < 1


@pytest.mark.parametrize("p_del,p_ins,p_sub", [(0.05, 0.0, 0.0), (0.0, 0.03, 0.0), (0.0, 0.0, 0.04)])
def test_empirical_rates_match(p_del, p_ins, p_sub):
    m = SynthesisModel("t", p_del, p_ins, p_sub, "test")
    rng = np.random.default_rng(3)
    s = _seqs(300, 200, seed=4)
    muts = [mutate(x, m, rng) for x in s]
    n = 300 * 200
    if p_sub:
        rate = sum(a != b for x, y in zip(s, muts) for a, b in zip(x, y)) / n
        assert abs(rate - p_sub * 1.0) < 0.006  # each substitution changes the base
    else:
        dl = sum(len(y) - len(x) for x, y in zip(s, muts)) / n
        assert abs(dl - (p_ins - p_del)) < 0.006


def test_context_independent_homopolymer_vs_alternating():
    """Documented property: no homopolymer dependence (none is supported by the literature)."""
    m = SynthesisModel("t", 0.05, 0.0, 0.0, "test")
    rng = np.random.default_rng(5)
    hp = "AAAA" * 50
    alt = "ACGT" * 50
    d_hp = np.mean([200 - len(mutate(hp, m, rng)) for _ in range(2000)])
    d_alt = np.mean([200 - len(mutate(alt, m, rng)) for _ in range(2000)])
    assert abs(d_hp - d_alt) < 0.5


def test_edit_distance_scale():
    m = PRESETS["photolitho_lietard2021"]
    rng = np.random.default_rng(6)
    s = _seqs(200, 150, seed=8)
    ed = np.mean([edlib.align(mutate(x, m, rng), x)["editDistance"] for x in s]) / 150
    assert 0.04 < ed < 0.08  # ~ 0.0465 + 0.0058 + 0.0098 minus alignment merging
