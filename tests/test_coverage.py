from __future__ import annotations

import numpy as np
import pytest

from dnastore.coverage import SENSITIVITY_GRID, CoverageModel, acc_factor, expected_copies, gc_factor, with_overrides


def _feat(n=2000, gc=0.5, acc=None):
    f = {"gc": np.full(n, gc)}
    if acc is not None:
        f["accessibility"] = np.full(n, acc)
    return f


def test_uniform_preset_is_flat():
    c = expected_copies(_feat(), SENSITIVITY_GRID["uniform"], np.random.default_rng(0))
    assert np.allclose(c, c[0])


def test_synthesis_only_reduces_to_lognormal():
    m = SENSITIVITY_GRID["twist_synth_only"]
    c = expected_copies(_feat(50_000), m, np.random.default_rng(1))
    logc = np.log(c / ((1 + m.e_bar) ** m.n_cycles * m.n0))
    assert abs(logc.std() - m.sigma_syn) < 0.01
    assert abs(np.mean(c / ((1 + m.e_bar) ** m.n_cycles * m.n0)) - 1.0) < 0.01  # mean-1 lognormal
    # identical draws to a plain lognormal with the same generator
    ref = m.n0 * np.random.default_rng(1).lognormal(-0.5 * m.sigma_syn**2, m.sigma_syn, 50_000) * (1 + m.e_bar) ** m.n_cycles
    assert np.allclose(c, ref)


@pytest.mark.parametrize("name", list(SENSITIVITY_GRID))
def test_presets_reproducible(name):
    f = _feat(500, acc=0.3)
    a = expected_copies(f, SENSITIVITY_GRID[name], np.random.default_rng(42))
    b = expected_copies(f, SENSITIVITY_GRID[name], np.random.default_rng(42))
    assert np.array_equal(a, b) and np.all(a >= 0)


def test_lower_efficiency_fewer_copies():
    base = CoverageModel(n0=100, n_cycles=20, e_bar=0.9)
    hi = expected_copies(_feat(10), base, np.random.default_rng(0))
    lo = expected_copies(_feat(10), with_overrides(base, e_bar=0.8), np.random.default_rng(0))
    assert np.all(lo < hi)
    # stochastic branching: lower efficiency -> lower mean
    sb = with_overrides(base, stochastic=True)
    m_hi = expected_copies(_feat(3000), sb, np.random.default_rng(0)).mean()
    m_lo = expected_copies(_feat(3000), with_overrides(sb, e_bar=0.8), np.random.default_rng(0)).mean()
    assert m_lo < m_hi


def test_stochastic_mean_matches_expectation():
    m = CoverageModel(n0=50, n_cycles=12, e_bar=0.9, stochastic=True)
    c = expected_copies(_feat(20_000), m, np.random.default_rng(3))
    assert abs(c.mean() / (50 * 1.9**12) - 1) < 0.02


def test_gc_factor_aird_monotone_flanks():
    g = gc_factor(np.array([0.05, 0.3, 0.5, 0.6, 0.7]), "aird2011_phusion_standard")
    assert g[1] == g[2] == 1.0
    assert g[3] < 1.0 and g[4] < g[3] and g[0] < 1.0
    assert np.isclose(g[4] ** 10, 0.01)
    assert np.all(gc_factor(np.array([0.2, 0.8]), "none") == 1.0)


def test_accessibility_penalty_monotone():
    f = acc_factor(np.array([0.0, 0.25, 0.5, 0.9]), penalty=0.2, threshold=0.5)
    assert np.allclose(f, [0.8, 0.9, 1.0, 1.0])
    assert acc_factor(None, 0.2, 0.5) is None
    m = SENSITIVITY_GRID["acc_assumed_strong"]
    m = with_overrides(m, sigma_syn=0.0, sd_rel=0.0, stochastic=False)
    open_ = expected_copies(_feat(10, acc=0.9), m, np.random.default_rng(0))
    closed = expected_copies(_feat(10, acc=0.0), m, np.random.default_rng(0))
    assert np.all(closed < open_)


def test_every_nondefault_preset_cites_sources():
    for name, m in SENSITIVITY_GRID.items():
        if name != "uniform":
            assert m.sources, name
        if m.acc_penalty > 0:
            assert any("ASSUMED" in v for v in m.sources.values())
