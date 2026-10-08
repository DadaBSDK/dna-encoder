"""Secondary-structure metrics: primer-site accessibility (primary), MFE (secondary), and a
cheap incremental proxy for use inside the steering beam search.

Definitions are in docs/notes/structure.md. Summary:

Strands
-------
The oligo as synthesised is the *top* strand ``T = F + body + rc(R)`` (length L). The
*bottom* strand is ``B = rc(T) = R + rc(body) + rc(F)``.
* Primer R anneals to T at its 3'-terminal ``r`` positions ``[L-r, L)`` (the site is rc(R)).
* Primer F anneals to B at its 3'-terminal ``f`` positions ``[L-f, L)`` (the site is rc(F)).

Accessibility
-------------
``P_u(i)`` = equilibrium probability that position i is unpaired, from the full
McCaskill partition function (ViennaRNA, DNA Mathews 2004 parameters, explicit
temperature and monovalent salt). ``site_access_X`` is the mean ``P_u`` over primer X's
site. The **3'-anchor**: primer and site are antiparallel, so the primer's 3'-terminal
nucleotide pairs with the 5'-most site position. ``anchor_access_X`` is therefore the
mean ``P_u`` over the first ``anchor_len`` positions of the site, i.e. ``[L-r, L-r+a)`` on
T for R and ``[L-f, L-f+a)`` on B for F. Polymerase extension starts there, so this is
the most functionally relevant region.

Proxy (incremental)
-------------------
See :func:`proxy_increments`. The full-sequence score is *defined* as the sum of
per-position increments, so an incremental beam implementation can be checked for exact
equality.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from functools import lru_cache
from typing import Any

import numpy as np

from .dnautil import revcomp

_DNA_LOADED = False


def _ensure_dna_params() -> None:
    """Load DNA (Mathews 2004) parameters once per process (safe in pool workers)."""
    global _DNA_LOADED
    if not _DNA_LOADED:
        import RNA

        RNA.params_load_DNA_Mathews2004()
        _DNA_LOADED = True


@dataclass
class StructureConfig:
    """Conditions and (calibrated) thresholds for structure metrics.

    ``salt_molar`` is monovalent cation concentration (mol/L), passed to ViennaRNA's
    ``md.salt``. ViennaRNA 2.7.2 has no Mg2+ model (see docs/notes/structure.md).
    Thresholds are ``None`` until calibrated; when set, they are lower bounds on
    accessibility at ``primary_temperature``.
    """

    temperatures: tuple[float, ...] = (60.0, 37.0)
    primary_temperature: float = 60.0
    salt_molar: float = 0.05
    anchor_len: int = 8
    min_anchor_access: float | None = None
    min_site_access: float | None = None

    @classmethod
    def from_dict(cls, d: dict | None) -> StructureConfig:
        d = dict(d or {})
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ValueError(f"unknown StructureConfig keys: {sorted(unknown)}")
        if "temperatures" in d:
            d["temperatures"] = tuple(float(t) for t in d["temperatures"])
        return cls(**d)


def make_md(temperature: float, salt: float) -> Any:
    """ViennaRNA model details with explicit temperature (deg C) and monovalent salt (M)."""
    import RNA

    _ensure_dna_params()
    md = RNA.md()
    md.temperature = float(temperature)
    md.salt = float(salt)
    return md


def mfe_at(seq: str, temperature: float, salt: float) -> float:
    """Minimum free energy (kcal/mol) at the given conditions."""
    import RNA

    fc = RNA.fold_compound(seq, make_md(temperature, salt))
    return float(fc.mfe()[1])


def unpaired_probs(seq: str, temperature: float, salt: float) -> np.ndarray:
    """Per-position probability of being unpaired, from the full partition function."""
    import RNA

    fc = RNA.fold_compound(seq, make_md(temperature, salt))
    _, e = fc.mfe()
    fc.exp_params_rescale(e)  # avoids overflow for long/structured sequences
    fc.pf()
    p = np.asarray(fc.bpp())[1:, 1:]
    p = p + p.T
    return np.clip(1.0 - p.sum(axis=1), 0.0, 1.0)


def accessibility(oligo: str, f_len: int, r_len: int, temperature: float, salt: float, anchor_len: int = 8) -> dict[str, float]:
    """Primer-site and 3'-anchor accessibility on both strands (see module docstring)."""
    L = len(oligo)
    pu_top = unpaired_probs(oligo, temperature, salt)
    pu_bot = unpaired_probs(revcomp(oligo), temperature, salt)
    a_r, a_f = min(anchor_len, r_len), min(anchor_len, f_len)
    out = {
        "site_access_R": float(pu_top[L - r_len :].mean()),
        "site_access_F": float(pu_bot[L - f_len :].mean()),
        "anchor_access_R": float(pu_top[L - r_len : L - r_len + a_r].mean()),
        "anchor_access_F": float(pu_bot[L - f_len : L - f_len + a_f].mean()),
    }
    out["min_site_access"] = min(out["site_access_R"], out["site_access_F"])
    out["min_anchor_access"] = min(out["anchor_access_R"], out["anchor_access_F"])
    return out


def structure_report(oligo: str, f_len: int, r_len: int, cfg: StructureConfig | None = None) -> dict[str, float]:
    """Accessibility and MFE (both strands) at every configured temperature.

    Keys are suffixed with ``_T{int(T)}``, e.g. ``min_anchor_access_T60`` and
    ``mfe_rc_T37``. If thresholds are set, ``access_ok`` reports whether the oligo passes
    at the primary temperature.
    """
    cfg = cfg or StructureConfig()
    rep: dict[str, float] = {}
    for t in cfg.temperatures:
        tag = f"_T{int(round(t))}"
        for k, v in accessibility(oligo, f_len, r_len, t, cfg.salt_molar, cfg.anchor_len).items():
            rep[k + tag] = v
        rep["mfe_fwd" + tag] = mfe_at(oligo, t, cfg.salt_molar)
        rep["mfe_rc" + tag] = mfe_at(revcomp(oligo), t, cfg.salt_molar)
    tag = f"_T{int(round(cfg.primary_temperature))}"
    if cfg.min_anchor_access is not None or cfg.min_site_access is not None:
        if "min_anchor_access" + tag not in rep:
            rep.update({k + tag: v for k, v in accessibility(oligo, f_len, r_len, cfg.primary_temperature,
                                                               cfg.salt_molar, cfg.anchor_len).items()})
        ok = True
        if cfg.min_anchor_access is not None:
            ok &= rep["min_anchor_access" + tag] >= cfg.min_anchor_access
        if cfg.min_site_access is not None:
            ok &= rep["min_site_access" + tag] >= cfg.min_site_access
        rep["access_ok"] = bool(ok)
    return rep


# --------------------------------------------------------------------------- proxy

# SantaLucia (1998) unified nearest-neighbour dG37 (kcal/mol) for Watson-Crick stacks,
# keyed by the 5'->3' dinucleotide on one strand.
_NN_DG37 = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}


@lru_cache(maxsize=None)
def stem_weight(kmer: str) -> float:
    """Stability weight of a k-mer duplex: ``-sum(NN dG37)`` over its k-1 stacks (> 0)."""
    return -sum(_NN_DG37[kmer[i : i + 2]] for i in range(len(kmer) - 1))


@dataclass(frozen=True)
class ProxyParams:
    """Parameters of the incremental structure proxy.

    ``ks``: k-mer lengths scanned. ``min_loop``: minimum hairpin loop between a k-mer and
    an earlier complementary k-mer. ``w_self`` / ``w_site``: weights of the self-
    complementarity term and the primer-site complementarity term.
    """

    ks: tuple[int, ...] = (5, 6)
    min_loop: int = 3
    w_self: float = 1.0
    w_site: float = 1.0


def site_kmer_sets(forward: str, reverse: str, ks: tuple[int, ...]) -> dict[int, frozenset[str]]:
    """k-mers of the top strand that are complementary to either primer site.

    R's site on T is rc(R) (the last r nt), so a T k-mer can pair with it iff it is a
    k-mer of R. F's site on B is rc(F); a T k-mer z appears on B as rc(z), which pairs
    with rc(F) iff z is a k-mer of rc(F). On T itself, z then also pairs with the F region.
    The term is therefore strand-symmetric.
    """
    F, R = forward, reverse
    rcF = revcomp(F)
    out = {}
    for k in ks:
        s = set()
        for src in (R, rcF):
            s.update(src[i : i + k] for i in range(len(src) - k + 1))
        out[k] = frozenset(s)
    return out


def proxy_increments(oligo: str, f_len: int, r_len: int, params: ProxyParams | None = None,
                     primers: tuple[str, str] | None = None) -> list[float]:
    """Per-position proxy increments. ``proxy_score == sum(proxy_increments)``.

    Incremental definition. When base ``i`` is appended, for each ``k`` in ``params.ks``
    with ``i >= k-1``, let ``z = oligo[i-k+1 : i+1]`` (the k-mer ending at i, starting at
    ``s = i-k+1``) and ``w = stem_weight(z)``. Then:

    * **self term**: ``w_self * w * C_k[rc(z)]``, where ``C_k`` counts k-mers that *end* at
      ``e <= s - 1 - min_loop`` (an earlier complementary k-mer with a loop of at least
      ``min_loop`` between them). Maintain ``C_k`` by adding, before scoring position i,
      the k-mer ending at ``e = i - k - min_loop`` (if ``e >= k-1``). That is O(1) per k
      per base.
    * **site term**: ``w_site * w * [z in S_k]``, where ``S_k`` is the fixed set of
      k-mers of ``R`` and of ``rc(F)`` (see :func:`site_kmer_sets`), precomputed once
      from the primers. Junction and primer-region k-mers are scored the same way; the
      constant primer-only part is identical across candidates and does not affect
      ranking.

    ``primers=(forward, reverse)`` fixes the site sets explicitly (needed when scoring a
    *prefix* during a beam search). By default they are read from the oligo itself:
    ``F = oligo[:f_len]``, ``R = rc(oligo[-r_len:])``.
    """
    p = params or ProxyParams()
    L = len(oligo)
    if primers is None:
        primers = (oligo[:f_len], revcomp(oligo[len(oligo) - r_len :]))
    sites = site_kmer_sets(primers[0], primers[1], p.ks)
    counts: dict[int, dict[str, int]] = {k: {} for k in p.ks}
    inc = [0.0] * L
    for i in range(L):
        total = 0.0
        for k in p.ks:
            e = i - k - p.min_loop
            if e >= k - 1:
                z_old = oligo[e - k + 1 : e + 1]
                counts[k][z_old] = counts[k].get(z_old, 0) + 1
            if i < k - 1:
                continue
            z = oligo[i - k + 1 : i + 1]
            w = stem_weight(z)
            c = counts[k].get(revcomp(z), 0)
            total += p.w_self * w * c + (p.w_site * w if z in sites[k] else 0.0)
        inc[i] = total
    return inc


def proxy_components(oligo: str, f_len: int, r_len: int, ks: tuple[int, ...] = (5, 6), min_loop: int = 3) -> dict[str, float]:
    """Self and site terms separately (for validation of combinations)."""
    self_ = sum(proxy_increments(oligo, f_len, r_len, ProxyParams(ks, min_loop, 1.0, 0.0)))
    site = sum(proxy_increments(oligo, f_len, r_len, ProxyParams(ks, min_loop, 0.0, 1.0)))
    return {"self": self_, "site": site}


def proxy_score(oligo: str, f_len: int, r_len: int, params: ProxyParams | None = None) -> float:
    """Full-sequence proxy score (higher = more predicted structure = worse)."""
    return float(sum(proxy_increments(oligo, f_len, r_len, params)))


def limited_span_mfe(seq: str, temperature: float, salt: float, max_bp_span: int) -> float:
    """Reference-only (not incremental): MFE with base pairs limited to ``max_bp_span``."""
    import RNA

    md = make_md(temperature, salt)
    md.max_bp_span = int(max_bp_span)
    return float(RNA.fold_compound(seq, md).mfe()[1])


def local_access(oligo: str, f_len: int, r_len: int, window: int, temperature: float, salt: float,
                 anchor_len: int = 8) -> dict[str, float]:
    """Tier-2 proxy: primer-site accessibility from partition functions of the terminal
    ``window`` nt only.

    * R site: the last ``window`` nt of the top strand T.
    * F site: the last ``window`` nt of the bottom strand B, which is ``rc(T[:window])``. It
      therefore depends only on the *first* ``window`` nt of the oligo, so a left-to-right
      beam can evaluate it exactly once it has placed ``window`` bases.

    Structure that pairs a site with bases farther than ``window`` nt away is ignored. That
    is the approximation; its rank agreement with :func:`accessibility` is reported in
    docs/notes/structure.md.
    """
    w = min(window, len(oligo))
    top = unpaired_probs(oligo[-w:], temperature, salt)
    bot = unpaired_probs(revcomp(oligo[:w]), temperature, salt)
    a_r, a_f = min(anchor_len, r_len), min(anchor_len, f_len)
    out = {
        "local_site_R": float(top[w - r_len :].mean()),
        "local_site_F": float(bot[w - f_len :].mean()),
        "local_anchor_R": float(top[w - r_len : w - r_len + a_r].mean()),
        "local_anchor_F": float(bot[w - f_len : w - f_len + a_f].mean()),
    }
    out["local_min_site"] = min(out["local_site_R"], out["local_site_F"])
    out["local_min_anchor"] = min(out["local_anchor_R"], out["local_anchor_F"])
    return out
