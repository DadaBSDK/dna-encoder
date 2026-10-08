"""NOVEL arm: steering-base constrained encoding (alphabets B and A).

Layout of the payload region
----------------------------
Position ``p`` (0-based within the payload region) is a **steering** position iff ``P`` is
finite and ``p % (P + 1) == P``, i.e. the pattern is ``[P info][1 steer][P info][1 steer]...``.
Steering bases carry no data. The decoder knows their positions and skips them.

* **Alphabet B (primary):** each info base is chosen from the 3 bases != the previous base
  (rotating ternary; the previous base may be a steering base). One steering choice
  therefore remaps the *whole* following info segment (DESIGN.md F5). Runs are at most 2 by
  construction. The frame is XOR-whitened with a keystream keyed by (global seed, index,
  seed) and carried as ``n_info`` trits by big-integer conversion.
* **Alphabet A (comparison arm):** info bases carry 2 bits each of the XOR-whitened frame.
  Info positions beyond ``4 * frame_bytes`` (when ``n_info % 4 != 0``) become extra steering
  positions.

Beam search
-----------
Left flank (primer F + index + seed), payload and right flank (revcomp R) are scored as
one string. At each steering position every beam state branches on the 4 bases, then the
deterministic info segment that follows is appended (for B it depends on the branch).
Incremental **hard** checks prune branches: run length, every completed GC window, exact
3'-anchored primer 8-mers in both orientations, and a feasibility bound on global GC.
Incremental **soft** costs rank them (weights are fixed a priori in the config): squared
window-GC deviation, k-mer confusability on both strands, and a self-complementarity stem
proxy. The final candidates get the full-oligo screen (``screening.screen_oligo``) and,
optionally, rescoring with the structure module. If no candidate meets every hard
constraint, the least-bad one is returned with ``info["violations"]`` set; the encoder
then rerolls the seed. Nothing is relaxed silently.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from ..dnautil import (
    BASES,
    base4_to_int,
    int_to_base4,
    int_to_trits,
    keystream_int,
    max_bits_for_trits,
    revcomp,
    trits_to_int,
)
from . import Codec, EncodedPayload, register

DOMAIN_A = 10
DOMAIN_B = 11
DOMAIN_C = 12
_BI = {b: i for i, b in enumerate(BASES)}
_GC = {"A": 0, "C": 1, "G": 1, "T": 0}


def steering_mask(payload_nt: int, P: int | None) -> list[bool]:
    """``True`` at steering positions of a payload region."""
    if not P:
        return [False] * payload_nt
    return [(p % (P + 1)) == P for p in range(payload_nt)]


# --------------------------------------------------------------------------- cost hooks
# Optional per-process tables, attached lazily so that codec objects stay cheap to pickle.

_CONF_CACHE: dict[tuple, tuple[np.ndarray, np.ndarray, int]] = {}


def _conf_tables(model: str, delta: float | None) -> tuple[np.ndarray, np.ndarray, int]:
    """(conf[code], rc_code[code], k) from :mod:`dnastore.kmer` (loaded once per process)."""
    key = (model, delta)
    if key not in _CONF_CACHE:
        from .. import kmer  # deferred: requires downloaded tables

        m = kmer.load_model(model)
        if m.units == "std":
            table = kmer.confusability_table(m, delta=delta) if delta is not None else kmer.confusability_table(m)
        else:
            table = kmer.confusability_table(m, c=delta) if delta is not None else kmer.confusability_table(m)
        _CONF_CACHE[key] = (np.asarray(table, dtype=np.float64), np.asarray(kmer.rc_code_table(m.k)), m.k)
    return _CONF_CACHE[key]


@dataclass
class SteeringOptions:
    """Encoder-side options (not stored in the header; recorded in the manifest).

    Constraint *values* (GC limits, run length, anchor length) are taken from the screening
    config passed to :meth:`Steering.bind`, so the beam and the screen cannot disagree.
    """

    beam_width: int = 16
    w_gc: float = 1.0
    w_conf: float = 1.0
    w_stem: float = 1.0  # overall weight of the tier-1 structure proxy
    stem_ks: tuple[int, ...] = (4, 5, 6)  # validated best (docs/notes/structure.md)
    stem_min_loop: int = 3
    stem_w_self: float = 1.0
    stem_w_site: float = 8.0
    conf_model: str = "r10.4.1_9mer"
    conf_delta: float | None = None  # None = model default
    conf_quantile: float = 0.75  # only confusability above this per-k-mer quantile is penalised
    rescore: str | None = "local_access"  # None | "local_access" (tier-2 proxy of the primary metric)
    rescore_top: int = 4
    access_window: int = 80
    access_temperature: float = 60.0
    access_salt: float = 0.05
    normalize: bool = True  # divide each soft weight by its term's IQR on whiten+RS oligos
    c_margin_sd: float = 4.0  # alphabet C capacity margin (see c_margin_bits)

    @classmethod
    def from_dict(cls, d: dict[str, Any] | None) -> SteeringOptions:
        d = dict(d or {})
        if "stem_ks" in d:
            d["stem_ks"] = tuple(d["stem_ks"])
        return cls(**d)


@dataclass
class _State:
    seq: str
    run: int
    gc_total: int
    gc_win: int  # GC count of the last ``gc_window_size`` bases
    cost: float
    ccode: int  # rolling code of the last conf k-mer
    stem: tuple  # per k: (scode, rcode, codes_by_end_position, counts)


_STEM_TABLES: dict[tuple, tuple] = {}


def _stem_tables(k: int, forward: str, reverse: str) -> tuple[list[float], bytearray]:
    """(stem_weight by code, site-membership by code) for the tier-1 proxy."""
    key = (k, forward, reverse)
    if key not in _STEM_TABLES:
        from ..structure import site_kmer_sets, stem_weight

        sites = site_kmer_sets(forward, reverse, (k,))[k]
        weights, member = [0.0] * 4**k, bytearray(4**k)
        for code in range(4**k):
            km = "".join(BASES[(code >> (2 * (k - 1 - j))) & 3] for j in range(k))
            weights[code] = stem_weight(km)
            member[code] = km in sites
        _STEM_TABLES[key] = (weights, member)
    return _STEM_TABLES[key]


class _Scorer:
    """Incremental hard checks and soft costs over ``left + payload + right``.

    Soft costs are additive over appended positions, so the final cost is an exact sum over
    the finished sequence. The left flank is appended with ``check=False``: its terms are
    constant across candidates and are not charged.

    * GC: ``w_gc * (g_w - 0.5)**2`` for every completed window.
    * Confusability: ``w_conf * (excess(x) + excess(rc(x)))`` per completed k-mer ``x``, with
      ``excess = max(0, conf - conf_floor)`` and ``conf_floor`` = the ``conf_quantile``
      quantile of the per-k-mer table. The per-oligo *mean* barely varies between oligos
      (docs/notes/confusability.md), so only the upper tail is targeted.
    * Structure, tier 1: ``w_stem * dnastore.structure.proxy_increments`` with
      ``ks = stem_ks``, ``w_self``, ``w_site``, ``min_loop`` (NN-weighted self-complementary
      k-mer pairs plus primer-site complementarity). The implementation is identical by
      construction and tested for equality.
    """

    def __init__(self, opts: SteeringOptions, scfg: Any, left: str, right: str, total_len: int, anchors: set[str],
                 primers: Any):
        self.o = opts
        self.c = scfg
        self.left, self.right, self.L = left, right, total_len
        self.anchors = anchors
        self.true_anchor_start = total_len - len(right)  # revcomp(R)[:a] at the R site is genuine
        self.conf = None
        if opts.w_conf:
            conf, rc_code, self.ck = _conf_tables(opts.conf_model, opts.conf_delta)
            floor = float(np.quantile(conf, opts.conf_quantile))
            self.conf = np.maximum(conf - floor, 0.0).astype(np.float64).tolist()
            self.rc_code = rc_code.tolist()
            self.cmask = (1 << (2 * self.ck)) - 1
        self.stem_spec = []
        if opts.w_stem:
            for k in opts.stem_ks:
                w, m = _stem_tables(k, primers.forward, primers.reverse)
                self.stem_spec.append((k, (1 << (2 * k)) - 1, 2 * (k - 1), w, m))

    def initial(self) -> _State:
        stem = tuple((0, 0, [], bytearray(4**k)) for k, *_ in self.stem_spec)
        st = _State("", 0, 0, 0, 0.0, 0, stem)
        return self.extend(st, self.left, check=False)

    def extend(self, st: _State, add: str, check: bool = True) -> _State | None:
        """Append bases; return the new state, or ``None`` if a hard constraint fails."""
        o, c = self.o, self.c
        seq, run, gct, gcw, cost, ccode = st.seq, st.run, st.gc_total, st.gc_win, st.cost, st.ccode
        w, a = c.gc_window_size, c.primer_anchor_len
        wlo, whi, maxhp = c.gc_window[0] * w, c.gc_window[1] * w, c.max_homopolymer
        conf, w_conf, w_gc = self.conf, o.w_conf, o.w_gc
        anchors, tas = self.anchors, self.true_anchor_start
        ml, w_self, w_site, w_stem = o.stem_min_loop, o.stem_w_self, o.stem_w_site, o.w_stem
        stem = [[sc, rc, list(codes), counts, False] for sc, rc, codes, counts in st.stem]
        comb = seq + add
        n = len(seq)
        last = seq[-1] if seq else ""
        for b in add:
            bi = _BI[b]
            run = run + 1 if b == last else 1
            last = b
            g = _GC[b]
            gct += g
            gcw += g
            if n >= w:
                gcw -= _GC[comb[n - w]]
            i = n  # 0-based position of this base
            n += 1
            if check:
                if run > maxhp:
                    return None
                if n >= w:
                    if gcw < wlo or gcw > whi:
                        return None
                    if w_gc:
                        cost += w_gc * (gcw / w - 0.5) ** 2
                if n >= a and n - a != tas and comb[n - a : n] in anchors:
                    return None
            if conf is not None:
                ccode = ((ccode << 2) | bi) & self.cmask
                if check and n >= self.ck:
                    cost += w_conf * (conf[ccode] + conf[self.rc_code[ccode]])
            for spec, stt in zip(self.stem_spec, stem):
                k, kmask, rshift, weights, member = spec
                stt[0] = ((stt[0] << 2) | bi) & kmask
                stt[1] = (stt[1] >> 2) | ((3 - bi) << rshift)
                codes = stt[2]
                codes.append(stt[0] if i >= k - 1 else -1)
                e = i - k - ml
                if e >= k - 1:
                    if not stt[4]:
                        stt[3] = bytearray(stt[3])
                        stt[4] = True
                    stt[3][codes[e]] += 1
                if check and i >= k - 1:
                    z = stt[0]
                    cost += w_stem * weights[z] * (w_self * stt[3][stt[1]] + (w_site if member[z] else 0.0))
        if check:
            remaining = self.L - n
            if gct + remaining < c.gc_global[0] * self.L or gct > c.gc_global[1] * self.L:
                return None
        return _State(comb, run, gct, gcw, cost, ccode, tuple((x[0], x[1], x[2], x[3]) for x in stem))


# --------------------------------------------------------------------------- normalisation

_SPREAD_CACHE: dict[str, dict[str, float]] = {}
CALIBRATION_N = 400


def _calibration_key(opts: SteeringOptions, scfg: Any, primers: Any, oligo_len: int, payload_nt: int, prefix_nt: int) -> str:
    import hashlib
    import json

    blob = json.dumps({
        "v": 1, "F": primers.forward, "R": primers.reverse, "L": oligo_len, "payload_nt": payload_nt,
        "prefix_nt": prefix_nt, "gc_window": scfg.gc_window_size, "conf": [opts.conf_model, opts.conf_delta,
        opts.conf_quantile], "stem": [list(opts.stem_ks), opts.stem_min_loop, opts.stem_w_self, opts.stem_w_site],
        "n": CALIBRATION_N,
    }, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def calibrate_spreads(opts: SteeringOptions, scfg: Any, primers: Any, left_len: int, payload_nt: int,
                      oligo_len: int) -> dict[str, Any]:
    """Spread (IQR, plus std for reference) of each soft-cost term's per-oligo total on
    unconstrained whitened quaternary oligos (the whiten+RS arm's sequence distribution).

    Each soft weight is divided by its term's spread (DESIGN.md Part 0c), so a weight of 1
    means "one baseline IQR of that term". Deterministic: fixed seed, ``CALIBRATION_N`` oligos.
    The term totals are computed by the beam's own :class:`_Scorer` with hard constraints
    disabled, so the calibration measures exactly what the beam optimises.
    """
    from ..screening import ScreeningConfig

    key = _calibration_key(opts, scfg, primers, oligo_len, payload_nt, left_len - len(primers.forward))
    if key in _SPREAD_CACHE:
        return _SPREAD_CACHE[key]
    relaxed = ScreeningConfig.from_dict({"gc_window_size": scfg.gc_window_size, "gc_window": [0.0, 1.0],
                                         "gc_global": [0.0, 1.0], "max_homopolymer": 10**6})
    rng = np.random.default_rng(20261002)
    terms = {"gc": dict(w_gc=1.0, w_conf=0.0, w_stem=0.0), "conf": dict(w_gc=0.0, w_conf=1.0, w_stem=0.0),
             "stem": dict(w_gc=0.0, w_conf=0.0, w_stem=1.0)}
    totals: dict[str, list[float]] = {t: [] for t in terms}
    scorers = {}
    for t, w in terms.items():
        o2 = SteeringOptions(**{**opts.__dict__, **w, "normalize": False})
        scorers[t] = o2
    for _ in range(CALIBRATION_N):
        left = primers.forward + "".join(BASES[x] for x in rng.integers(0, 4, left_len - len(primers.forward)))
        payload = "".join(BASES[x] for x in rng.integers(0, 4, payload_nt))
        total = len(left) + payload_nt + len(primers.tail)
        for t, o2 in scorers.items():
            sc = _Scorer(o2, relaxed, left, primers.tail, total, set(), primers)
            st = sc.extend(sc.initial(), payload)
            st = sc.extend(st, primers.tail)
            totals[t].append(st.cost)
    out: dict[str, Any] = {"n": CALIBRATION_N, "key": key}
    for t, v in totals.items():
        arr = np.asarray(v)
        q25, q75 = np.percentile(arr, [25, 75])
        iqr = float(q75 - q25)
        out[t] = {"iqr": iqr, "std": float(arr.std()), "median": float(np.median(arr))}
        out[t]["spread"] = iqr if iqr > 0 else (float(arr.std()) or 1.0)
    _SPREAD_CACHE[key] = out
    return out


def c_margin_bits(n_info: int, sd: float = 4.0) -> int:
    """Alphabet C safety margin: expected forced (radix-3) positions n/16 plus ``sd`` standard
    deviations, each costing 2 - log2(3) bits."""
    lam = n_info / 16.0
    return int(math.ceil((2 - math.log2(3)) * (lam + sd * math.sqrt(lam))))


@register
class Steering(Codec):
    """Steering-base codec. ``alphabet``: "A" (quaternary + whitening + seeds), "B" (rotating
    ternary), "C" (run-length-limited quaternary, max run = screening max_homopolymer)."""

    codec_id = 4
    name = "steering"

    def __init__(self, P: int | None = 8, alphabet: str = "A", **options: Any) -> None:
        if alphabet not in ("A", "B", "C"):
            raise ValueError("alphabet must be 'A', 'B' or 'C'")
        self.P = int(P) if P else None
        if self.P is not None and self.P < 1:
            raise ValueError("P must be nonnegative (0 means no steering)")
        self.alphabet = alphabet
        self.opts = SteeringOptions.from_dict(options)
        if self.opts.beam_width < 1 or self.opts.rescore_top < 1:
            raise ValueError("beam_width and rescore_top must be positive")
        self.primers = None
        self.scfg = None
        self.spreads: dict[str, Any] | None = None

    def bind(self, primers: Any, screening: dict[str, Any]) -> None:
        from ..screening import ScreeningConfig

        self.primers = primers
        self.scfg = ScreeningConfig.from_dict(screening)
        a = self.scfg.primer_anchor_len
        rc_f, rc_r = revcomp(primers.forward), revcomp(primers.reverse)
        # Exact 3'-anchored primer 8-mers in the top strand, as in screening.primer_hits:
        # F[-a:], R[-a:], revcomp(F)[:a], revcomp(R)[:a].
        self.anchors = {primers.forward[-a:], primers.reverse[-a:], rc_f[:a], rc_r[:a]}
        if self.scfg.primer_anchor_mismatch != 0:
            raise NotImplementedError("beam anchor check supports exact (0-mismatch) anchors only")

    # ------------------------------------------------------------------ geometry
    def _layout(self, payload_nt: int) -> tuple[list[bool], int]:
        mask = steering_mask(payload_nt, self.P)
        return mask, sum(1 for s in mask if not s)

    def frame_bits(self, payload_nt: int) -> int:
        n_info = self._layout(payload_nt)[1]
        if self.alphabet == "A":
            return 2 * n_info
        if self.alphabet == "B":
            return max_bits_for_trits(n_info)
        return 2 * n_info - c_margin_bits(n_info, self.opts.c_margin_sd)

    def _mask_bits(self, v: int, nbits: int, index: int, seed: int, global_seed: int) -> int:
        dom = {"A": DOMAIN_A, "B": DOMAIN_B, "C": DOMAIN_C}[self.alphabet]
        return v ^ keystream_int(nbits, global_seed, dom, index, seed)

    def _effective_opts(self, left_len: int, payload_nt: int) -> SteeringOptions:
        o = self.opts
        if not o.normalize:
            return o
        if self.spreads is None:
            total = left_len + payload_nt + len(self.primers.tail)
            self.spreads = calibrate_spreads(o, self.scfg, self.primers, left_len, payload_nt, total)
        sp = self.spreads
        return SteeringOptions(**{**o.__dict__, "w_gc": o.w_gc / sp["gc"]["spread"],
                                  "w_conf": o.w_conf / sp["conf"]["spread"],
                                  "w_stem": o.w_stem / sp["stem"]["spread"], "normalize": False})

    # ------------------------------------------------------------------ encode
    def encode_payload(self, frame: int, index: int, seed: int, left: str, right: str, payload_nt: int,
                       global_seed: int) -> EncodedPayload:
        if self.scfg is None:
            raise RuntimeError("Steering.bind(primers, screening) must be called before encoding")
        nbits = self.frame_bits(payload_nt)
        if frame >> nbits:
            raise ValueError("frame larger than capacity")
        mask, n_info = self._layout(payload_nt)
        v = self._mask_bits(frame, nbits, index, seed, global_seed)
        if self.alphabet == "A":
            symbols: list[int] | None = [_BI[b] for b in int_to_base4(v, n_info)]
        elif self.alphabet == "B":
            symbols = int_to_trits(v, n_info)
        else:
            symbols = None  # C: digits are extracted per beam state (radix depends on the run state)
            if len(left) >= 2 and left[-1] == left[-2]:
                raise RuntimeError("alphabet C assumes the prefix ends in a run of length 1")
        o = self._effective_opts(len(left), payload_nt)
        total = len(left) + payload_nt + len(right)
        sc = _Scorer(o, self.scfg, left, right, total, self.anchors, self.primers)
        maxrun = self.scfg.max_homopolymer

        # beam entries: (state, remaining value for C)
        beam: list[tuple[_State, int]] = [(sc.initial(), v)]
        sym_i, p, expansions = 0, 0, 0
        info_left = n_info
        while p < payload_nt:
            new: list[tuple[_State, int]] = []
            if mask[p]:
                for st, rem in beam:
                    for b in BASES:
                        nxt = sc.extend(st, b)
                        expansions += 1
                        if nxt is not None:
                            new.append((nxt, rem))
                p += 1
            else:
                q = p
                while q < payload_nt and not mask[q]:
                    q += 1
                seglen = q - p
                info_left -= seglen
                for st, rem in beam:
                    if symbols is not None:
                        seg, rem2 = self._render(symbols[sym_i : sym_i + seglen], st.seq[-1]), rem
                    else:
                        seg, rem2 = self._render_c(rem, seglen, st.seq[-1], st.run, maxrun)
                        if rem2 >> (2 * info_left):  # cannot fit the rest even at radix 4 everywhere
                            continue
                    nxt = sc.extend(st, seg)
                    expansions += 1
                    if nxt is not None:
                        new.append((nxt, rem2))
                sym_i += seglen
                p = q
            if not new:
                return self._fallback(v, symbols, mask, left, payload_nt, maxrun, "beam_exhausted", expansions)
            new.sort(key=lambda e: e[0].cost)
            seen, beam = set(), []
            for e in new:
                if e[0].seq not in seen:
                    seen.add(e[0].seq)
                    beam.append(e)
                if len(beam) >= o.beam_width:
                    break
        if symbols is None:
            beam = [e for e in beam if e[1] == 0]
            if not beam:
                return self._fallback(v, symbols, mask, left, payload_nt, maxrun, "c_overflow", expansions)
        return self._finalize([st for st, _ in beam], sc, left, right, payload_nt, expansions)

    def _render(self, syms: list[int], prev: str) -> str:
        if self.alphabet == "A":
            return "".join(BASES[s] for s in syms)
        out, pi = [], _BI[prev]
        for t in syms:
            pi = (pi + 1 + t) % 4
            out.append(BASES[pi])
        return "".join(out)

    @staticmethod
    def _render_c(v: int, n: int, prev: str, run: int, maxrun: int) -> tuple[str, int]:
        """Alphabet C: emit ``n`` info bases from value ``v`` (LSB-first mixed radix).
        Radix 3 (bases != prev, in ACGT order) when the current run equals ``maxrun``, else radix 4."""
        out, last = [], prev
        for _ in range(n):
            if run >= maxrun:
                v, d = divmod(v, 3)
                b = [x for x in BASES if x != last][d]
            else:
                v, d = divmod(v, 4)
                b = BASES[d]
            run = run + 1 if b == last else 1
            last = b
            out.append(b)
        return "".join(out), v

    def _fallback(self, v, symbols, mask, left, payload_nt, maxrun, reason, expansions) -> EncodedPayload:
        """Best-effort sequence (steering base = a base differing from its predecessor)."""
        seq, prev, run, si, rem = [], left[-1], 1, 0, v
        for p in range(payload_nt):
            if mask[p]:
                b = next(x for x in BASES if x != prev)
            elif symbols is not None:
                b = self._render([symbols[si]], prev)
                si += 1
            else:
                b, rem = self._render_c(rem, 1, prev, run, maxrun)
            run = run + 1 if b == prev else 1
            seq.append(b)
            prev = b
        viol = [reason] + (["c_overflow"] if symbols is None and rem and reason != "c_overflow" else [])
        return EncodedPayload("".join(seq), {"violations": viol, "expansions": expansions})

    def _finalize(self, beam, sc: _Scorer, left, right, payload_nt, expansions) -> EncodedPayload:
        """Full-oligo checks on the final beam, then tier-2 structure rescoring.

        Candidates are ordered by (has violations, beam cost). With ``rescore="local_access"``
        the ``rescore_top`` best violation-free candidates are re-ranked by the tier-2 proxy
        of the primary structure metric: min 3'-anchor accessibility from windowed partition
        functions (``structure.local_access``, window ``access_window``). Ties are broken by
        beam cost.
        """
        from ..screening import screen_oligo

        o = self.opts
        scored = []
        for st in beam:
            payload = st.seq[len(left) :]
            full = sc.extend(st, right)  # junction windows / runs / anchors / global GC / proxy
            if full is None:
                scored.append((True, math.inf, payload, ["final_check"], None))
                continue
            scr = screen_oligo(full.seq, self.primers, self.scfg, compute_mfe=False)
            scored.append((bool(scr["violations"]), full.cost, payload, scr["violations"], full.seq))
        scored.sort(key=lambda x: (x[0], x[1]))
        info: dict[str, Any] = {"expansions": expansions}
        if o.rescore == "local_access" and not scored[0][0]:
            from ..structure import local_access

            lf, lr = len(self.primers.forward), len(self.primers.reverse)
            top = [x for x in scored if not x[0]][: o.rescore_top]
            rescored = []
            for x in top:
                la = local_access(x[4], lf, lr, o.access_window, o.access_temperature, o.access_salt,
                                  self.scfg.primer_anchor_len)
                rescored.append((-la["local_min_anchor"], x[1], x, la["local_min_anchor"]))
            rescored.sort(key=lambda r: (r[0], r[1]))
            chosen, access = rescored[0][2], rescored[0][3]
            info["local_min_anchor"] = access
        else:
            chosen = scored[0]
        _, cost, payload, viol, _ = chosen
        return EncodedPayload(payload, {"violations": viol, "cost": cost, **info})

    # ------------------------------------------------------------------ decode
    def decode_payload(self, seq: str, index: int, seed: int, prev: str, global_seed: int) -> int | None:
        mask, n_info = self._layout(len(seq))
        nbits = self.frame_bits(len(seq))
        pi, run = _BI.get(prev), 1
        if pi is None:
            return None
        maxrun = self.scfg.max_homopolymer if self.scfg is not None else 3
        if self.alphabet == "C":
            v, mult = 0, 1
        else:
            syms: list[int] = []
        for p, b in enumerate(seq):
            q = _BI.get(b)
            if q is None:
                return None
            if not mask[p]:
                if self.alphabet == "B":
                    t = (q - pi - 1) % 4
                    if t == 3:
                        return None
                    syms.append(t)
                elif self.alphabet == "A":
                    syms.append(q)
                else:
                    if run >= maxrun:
                        if q == pi:
                            return None  # run would exceed the limit: impossible codeword
                        d = [i for i in range(4) if i != pi].index(q)
                        v += d * mult
                        mult *= 3
                    else:
                        v += q * mult
                        mult *= 4
            run = run + 1 if q == pi else 1
            pi = q
        if self.alphabet == "B":
            v = trits_to_int(syms)
        elif self.alphabet == "A":
            v = base4_to_int("".join(BASES[x] for x in syms))
        if v >> nbits:
            return None
        return self._mask_bits(v, nbits, index, seed, global_seed)

    def params(self) -> dict[str, Any]:
        p = {"P": self.P or 0, "alphabet": self.alphabet}
        if self.alphabet == "C":
            p["c_margin_sd"] = self.opts.c_margin_sd
            p["maxrun"] = self.scfg.max_homopolymer if self.scfg is not None else 3
        return p

    @classmethod
    def from_params(cls, params: dict[str, Any]) -> "Steering":
        params = dict(params)
        maxrun = params.pop("maxrun", None)
        obj = cls(**params)
        if maxrun is not None:
            from ..screening import ScreeningConfig

            obj.scfg = ScreeningConfig.from_dict({"max_homopolymer": maxrun})
        return obj
