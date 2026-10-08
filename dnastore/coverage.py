"""Optional sequence-dependent coverage model (synthesis + PCR amplification yield).

Purpose
-------
squigulator + Dorado model *sequencing* only. Arms that allow constraint violations
(``whiten+RS``) could pay a cost upstream, in PCR yield, which a sequencing simulator
never sees. This module gives each oligo a copy number that may depend on its GC
content and its primer-site accessibility. It is run **only as a sensitivity sweep**
over literature-derived presets (``SENSITIVITY_GRID``). See ``docs/notes/pcr_bias.md``
for every source, the exact location of each number, and the limitations.

Model
-----
For oligo ``i``:

1. Synthesis copies: ``s_i = n0 * L_i`` with ``L_i ~ LogNormal(-sigma_syn**2 / 2, sigma_syn)``
   (mean 1). Lognormal synthesis coverage with shape ``sigma_syn`` follows Gimpel et al. 2023.
2. Relative per-cycle efficiency (Gimpel et al. 2023 definition):
   ``r_i = (1 + e_i) / (1 + e_bar)``, with
   ``r_i = (1 + N(0, sd_rel)) * g_gc(GC_i) * g_acc(a_i) * tail_i``.
   * ``sd_rel``: spread of relative efficiencies (Gimpel 2023; Heckel 2019 illustrative).
   * ``g_gc``: 1 for presets backed by DNA-storage data (Chen 2020; Gimpel 2025: GC not
     predictive). The Aird 2011 curve is an out-of-domain pessimistic bound.
   * ``g_acc``: **assumed functional form, not literature-derived**. Only its magnitude
     range is anchored to Gimpel et al. 2025 effect sizes (4.8% motif effect, ~20% for the
     worst ~2% of sequences).
   * ``tail_i``: Gimpel 2025's poorly amplifying tail (~2% of sequences at relative
     efficiency as low as 0.8).
3. Absolute efficiency ``e_i = clip(r_i * (1 + e_bar) - 1, 0, 1)``. Amplify for ``n_cycles``:
   deterministically (``s_i * (1 + e_i)**c``), or stochastically with Chen et al. 2020's
   branching process ``n_{j+1} = n_j + Binomial(n_j, e_i)``. The exact binomial is used while
   ``n_j < switch_count`` and the expectation afterwards, which loses negligible noise.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np


@dataclass(frozen=True)
class CoverageModel:
    """Parameters of the synthesis + PCR coverage model. ``sources`` tags every non-default value."""

    name: str = "custom"
    sigma_syn: float = 0.0  # lognormal shape of synthesis coverage
    n0: float = 100.0  # mean molecules per oligo entering PCR (sampling choice, not a lit. value)
    n_cycles: int = 15
    e_bar: float = 0.95  # mean absolute per-cycle efficiency
    sd_rel: float = 0.0  # sd of relative efficiency (1+e_i)/(1+e_bar)
    gc_curve: str = "none"  # "none" | "aird2011_phusion_standard"
    acc_penalty: float = 0.0  # max fractional drop of relative efficiency at zero accessibility (ASSUMED form)
    acc_threshold: float = 0.5  # accessibility below which the penalty starts (ASSUMED)
    tail_frac: float = 0.0  # fraction of poorly amplifying oligos
    tail_rel_eff: float = 1.0  # their relative efficiency
    stochastic: bool = False
    switch_count: float = 1e4
    sources: dict[str, str] = field(default_factory=dict)


# Aird et al. 2011 (Genome Biol 12:R18), standard Phusion HF protocol, 10 cycles, Fig. 1e text:
# GC > 65% -> ~1/100 of mid-GC, GC < 12% -> ~1/10, plateau (>= 0.7) from 11% to 56% GC.
# Converted to per-cycle relative factors: 0.01**(1/10) and 0.1**(1/10). The plateau is set
# to 1 and the 56-65% and 11-12% flanks are log-linearly interpolated (ASSUMED shape between
# the reported anchors).
_AIRD_GC = np.array([0.0, 0.11, 0.12, 0.56, 0.65, 1.0])
_AIRD_LOGF = np.log(np.array([0.1 ** 0.1, 0.1 ** 0.1, 1.0, 1.0, 0.01 ** 0.1, 0.01 ** 0.1]))


def gc_factor(gc: np.ndarray, curve: str) -> np.ndarray:
    """Per-cycle relative efficiency multiplier as a function of oligo GC fraction."""
    gc = np.asarray(gc, dtype=float)
    if curve == "none":
        return np.ones_like(gc)
    if curve == "aird2011_phusion_standard":
        # The plateau anchors at 0.11 and 0.12 are within a percent; interpolation is fine.
        return np.exp(np.interp(gc, _AIRD_GC, _AIRD_LOGF))
    raise ValueError(f"unknown gc_curve {curve!r}")


def acc_factor(acc: np.ndarray | None, penalty: float, threshold: float) -> np.ndarray | None:
    """ASSUMED form: linear drop of relative efficiency below ``threshold`` accessibility,
    reaching ``1 - penalty`` at zero accessibility."""
    if acc is None or penalty == 0.0:
        return None
    acc = np.clip(np.asarray(acc, dtype=float), 0.0, 1.0)
    deficit = np.clip((threshold - acc) / threshold, 0.0, 1.0)
    return 1.0 - penalty * deficit


def relative_efficiency(features: dict[str, Any], model: CoverageModel, rng: np.random.Generator) -> np.ndarray:
    """Relative per-cycle efficiency ``r_i`` for every oligo."""
    gc = np.asarray(features["gc"], dtype=float)
    n = gc.size
    r = 1.0 + (rng.normal(0.0, model.sd_rel, n) if model.sd_rel > 0 else np.zeros(n))
    r *= gc_factor(gc, model.gc_curve)
    fa = acc_factor(features.get("accessibility"), model.acc_penalty, model.acc_threshold)
    if fa is not None:
        if fa.size != n:
            raise ValueError("accessibility array length mismatch")
        r *= fa
    if model.tail_frac > 0:
        r = np.where(rng.random(n) < model.tail_frac, r * model.tail_rel_eff, r)
    return r


def expected_copies(features: dict[str, Any], model: CoverageModel, rng: np.random.Generator) -> np.ndarray:
    """Molecule count of each oligo after synthesis and ``n_cycles`` of PCR.

    ``features`` must contain ``"gc"`` (array of GC fractions). It may contain
    ``"accessibility"`` (array in [0, 1], e.g. the minimum unpaired probability over the
    primer sites, computed by the caller). Normalise the result to obtain sampling weights.
    """
    gc = np.asarray(features["gc"], dtype=float)
    n = gc.size
    if model.sigma_syn > 0:
        s = model.n0 * rng.lognormal(-0.5 * model.sigma_syn**2, model.sigma_syn, n)
    else:
        s = np.full(n, float(model.n0))
    r = relative_efficiency(features, model, rng)
    e = np.clip(r * (1.0 + model.e_bar) - 1.0, 0.0, 1.0)
    if not model.stochastic:
        return s * (1.0 + e) ** model.n_cycles
    cnt = rng.poisson(s).astype(float)
    for _ in range(model.n_cycles):
        small = cnt < model.switch_count
        grow = np.empty_like(cnt)
        grow[small] = rng.binomial(cnt[small].astype(np.int64), e[small])
        grow[~small] = cnt[~small] * e[~small]
        cnt = cnt + grow
    return cnt


def with_overrides(model: CoverageModel, **kw: Any) -> CoverageModel:
    return replace(model, **kw)


# --------------------------------------------------------------------------- presets
_G23 = "Gimpel et al. 2023, doi:10.1038/s41467-023-41729-1"
_G25 = "Gimpel et al. 2025, doi:10.1038/s41467-025-64221-4"
_C20 = "Chen et al. 2020, doi:10.1038/s41467-020-16958-3"
_H19 = "Heckel et al. 2019, doi:10.1038/s41598-019-45832-6 (arXiv:1803.03322v1 Sec. 4)"
_A11 = "Aird et al. 2011, doi:10.1186/gb-2011-12-2-r18"
_ASSUMED = "ASSUMED - not literature-derived (magnitude anchored to Gimpel 2025 effect sizes)"

_BASE = CoverageModel(
    name="base",
    n_cycles=15,
    e_bar=0.95,
    sources={"n_cycles": f"{_G23} (15 cycles per round, six rounds)", "e_bar": f"{_C20} Fig. 5b model P=0.95 (origin not stated)"},
)

SENSITIVITY_GRID: dict[str, CoverageModel] = {
    # No bias at all: uniform copies (control).
    "uniform": replace(_BASE, name="uniform", sources={}),
    # Synthesis only, measured on Twist (material deposition), unconstrained pool.
    "twist_synth_only": replace(_BASE, name="twist_synth_only", sigma_syn=0.27,
                                sources={**_BASE.sources, "sigma_syn": f"{_G23} Fig. 2b (Twist, unconstrained)"}),
    # Twist synthesis + measured PCR efficiency spread; GC effect absent (Chen 2020, Gimpel 2025).
    "twist_pcr_measured": replace(_BASE, name="twist_pcr_measured", sigma_syn=0.27, sd_rel=0.0051, stochastic=True,
                                  sources={**_BASE.sources, "sigma_syn": f"{_G23} Fig. 2b",
                                           "sd_rel": f"{_G23} Fig. 3d (0.0051 unconstrained)",
                                           "gc_curve": f"none: {_C20} Fig. 4 (slope <0.01); {_G25} (GC not predictive)"}),
    # Upper end of the PCR efficiency spread across literature datasets.
    "pcr_spread_lit_upper": replace(_BASE, name="pcr_spread_lit_upper", sigma_syn=0.27, sd_rel=0.012, stochastic=True,
                                    sources={**_BASE.sources, "sd_rel": f"{_G23} Fig. 3d caption (range 0.0058-0.012)"}),
    # Electrochemical synthesis (CustomArray), unconstrained pool: much wider synthesis skew.
    "customarray_pcr": replace(_BASE, name="customarray_pcr", sigma_syn=1.30, sd_rel=0.0051, stochastic=True,
                               sources={**_BASE.sources, "sigma_syn": f"{_G23} Fig. 2b (CustomArray, unconstrained)"}),
    # Poorly amplifying tail (~2% of sequences at relative efficiency down to 0.8).
    "poor_tail": replace(_BASE, name="poor_tail", sigma_syn=0.27, sd_rel=0.0051, tail_frac=0.02, tail_rel_eff=0.8,
                         stochastic=True,
                         sources={**_BASE.sources, "tail": f"{_G25} Fig. 2c (~2% of pool, as low as 80%; lower bound used)"}),
    # Heckel's illustrative PCR model (a modelling assumption in that paper, citing literature).
    "heckel_illustrative": replace(_BASE, name="heckel_illustrative", e_bar=0.85, sd_rel=0.07 / 1.85, n_cycles=22,
                                   sources={"e_bar,sd_rel,n_cycles": f"{_H19}: factor ~ N(1.85, 0.07), 22 cycles (simulation assumption)"}),
    # Out-of-domain pessimistic GC bound (genomic libraries, Phusion standard protocol).
    "aird_gc_pessimistic": replace(_BASE, name="aird_gc_pessimistic", sigma_syn=0.27, sd_rel=0.0051,
                                   gc_curve="aird2011_phusion_standard", stochastic=True,
                                   sources={**_BASE.sources, "gc_curve": f"{_A11} Fig. 1e text (10 cycles, Phusion HF standard)"}),
    # Accessibility knob: ASSUMED form, two magnitudes anchored to Gimpel 2025 effect sizes.
    "acc_assumed_mild": replace(_BASE, name="acc_assumed_mild", sigma_syn=0.27, sd_rel=0.0051, acc_penalty=0.048,
                                stochastic=True,
                                sources={**_BASE.sources, "acc_penalty": f"{_ASSUMED}: 4.8% = TCGTGT insertion effect, {_G25} Fig. 6c"}),
    "acc_assumed_strong": replace(_BASE, name="acc_assumed_strong", sigma_syn=0.27, sd_rel=0.0051, acc_penalty=0.20,
                                  stochastic=True,
                                  sources={**_BASE.sources, "acc_penalty": f"{_ASSUMED}: 20% = worst ~2% sequences, {_G25} Fig. 2c"}),
}
