"""Synthesis error model: CONTEXT-FREE presets only (see docs/notes/synthesis_errors.md).

A literature search (Heckel et al. 2019; Gimpel et al. 2023, 2024; Lietard et al. 2021;
Yeom et al. 2023; Sabary et al. SOLQC; Wang et al. 2023) found **no published,
synthesis-specific, quantitative dependence of error rates on homopolymer run length or
GC content** in the regime this project uses (max run <= 3, GC 40-60 %). Only aggregate
per-nucleotide rates (and some position or base-identity effects) are published. This
module therefore applies **sequence-context-independent** per-nucleotide deletion,
insertion and substitution rates. It cannot penalise constrained vs. unconstrained
sequences differently, and conclusions about that comparison are scoped to sequencing +
amplification.

Every preset parameter carries a source tag. ``None`` rates mean "not reported" and are
treated as 0 with that fact recorded.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_BASES = "ACGT"


@dataclass(frozen=True)
class SynthesisModel:
    """Per-nucleotide, position- and context-independent synthesis error rates."""

    name: str
    p_del: float
    p_ins: float
    p_sub: float
    source: str
    notes: str = ""

    def __post_init__(self) -> None:
        for v in (self.p_del, self.p_ins, self.p_sub):
            if not 0.0 <= v < 1.0:
                raise ValueError("rates must lie in [0, 1)")


ZERO = SynthesisModel("zero", 0.0, 0.0, 0.0, "identity (no synthesis errors)")

PRESETS: dict[str, SynthesisModel] = {
    "zero": ZERO,
    "twist_gimpel2023": SynthesisModel(
        "twist_gimpel2023", 0.58e-3, 0.0, 0.0,
        "Gimpel et al. 2023, Nat Commun 14:6026, doi:10.1038/s41467-023-41729-1, Fig. 2a "
        "(material deposition, Twist): deletions 0.58 +/- 0.15 x 10^-3 nt^-1",
        "substitutions reported as negligible; synthesis insertions not reported (set 0); "
        "no positional dependence observed for material deposition",
    ),
    "electrochemical_gimpel2023": SynthesisModel(
        "electrochemical_gimpel2023", 13.5e-3, 0.0, 0.0,
        "Gimpel et al. 2023, doi:10.1038/s41467-023-41729-1, Fig. 2a (electrochemical, "
        "CustomArray/GenScript): deletions 13.5 +/- 2.0 x 10^-3 nt^-1",
        "substitutions negligible; insertions not reported (set 0). The paper reports a strong "
        "increase of deletions towards the 5' end and deletion clustering (10-14 % of deletions "
        "in runs); neither is modelled here (uniform rate = mean only)",
    ),
    "photolitho_gimpel2024": SynthesisModel(
        "photolitho_gimpel2024", 0.082, 0.016, 0.025,
        "Gimpel et al. 2024, bioRxiv doi:10.1101/2024.07.04.602085 (photolithographic "
        "synthesis): 0.082 deletions, 0.016 insertions, 0.025 substitutions per nt",
        "preprint; consecutive-error clustering (16 % of subs, 14 % of dels) not modelled",
    ),
    "photolitho_lietard2021": SynthesisModel(
        "photolitho_lietard2021", 0.0465, 0.0058, 0.0098,
        "Lietard et al. 2021, Nucleic Acids Res 49:6687, doi:10.1093/nar/gkab505, Table 1 "
        "(2SZ library): deletions 4.65 %, insertions 0.58 %, substitutions 0.98 % per bp",
        "G->T dominates substitutions (0.32 % per bp); base-specific substitution spectrum and "
        "the every-fifth-nucleotide deletion pattern (Fig. 6) are not modelled",
    ),
}


def mutate(seq: str, model: SynthesisModel, rng: np.random.Generator) -> str:
    """One synthesised molecule of ``seq``.

    Each template base is independently deleted (``p_del``), or else substituted
    (``p_sub``; uniform over the other three bases). After each template position an
    extra random base is inserted with probability ``p_ins``. Context-independent by
    construction.
    """
    if model.p_del == model.p_ins == model.p_sub == 0.0:
        return seq
    n = len(seq)
    u_del = rng.random(n)
    u_sub = rng.random(n)
    u_ins = rng.random(n)
    sub_off = rng.integers(1, 4, size=n)
    ins_base = rng.integers(0, 4, size=n)
    out: list[str] = []
    for i, b in enumerate(seq):
        if u_del[i] >= model.p_del:
            if u_sub[i] < model.p_sub:
                out.append(_BASES[(_BASES.index(b) + sub_off[i]) % 4])
            else:
                out.append(b)
        if u_ins[i] < model.p_ins:
            out.append(_BASES[ins_base[i]])
    return "".join(out)


def apply_synthesis_errors(oligos: list[str], model: SynthesisModel, rng: np.random.Generator,
                           copies: int | list[int] = 1) -> list[list[str]]:
    """Synthesised molecules for every oligo: ``copies`` molecules each (an int, or one count per oligo)."""
    counts = [copies] * len(oligos) if isinstance(copies, int) else list(copies)
    if len(counts) != len(oligos):
        raise ValueError("copies must match oligos")
    return [[mutate(s, model, rng) for _ in range(c)] for s, c in zip(oligos, counts)]
