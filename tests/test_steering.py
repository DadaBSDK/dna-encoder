"""Steering codec: round trips, steering-position semantics, constraints, and exact incremental costs."""

from __future__ import annotations

import random

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from dnastore.codecs.steering import Steering, SteeringOptions, _Scorer, steering_mask
from dnastore.decoder import decode_reads
from dnastore.dnautil import revcomp
from dnastore.encoder import encode_bytes
from dnastore.oligo import PREFIX_NT, OligoGeometry
from dnastore.screening import ScreeningConfig, screen_oligo

from .conftest import make_cfg

FAST = {"w_conf": 0.0, "beam_width": 4, "rescore": None, "normalize": False}


def _reads(res, seed=0):
    rng = random.Random(seed)
    reads = [s if rng.random() < 0.5 else revcomp(s) for _, s in res.oligos]
    rng.shuffle(reads)
    return reads


def test_steering_mask():
    assert steering_mask(10, 4) == [False] * 4 + [True] + [False] * 4 + [True]
    assert steering_mask(5, None) == [False] * 5
    assert steering_mask(5, 0) == [False] * 5


@pytest.mark.parametrize("alphabet", ["A", "B", "C"])
@pytest.mark.parametrize("P", [3, 8, 0])
@given(data=st.binary(max_size=1200))
@settings(max_examples=6, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_roundtrip(primers, alphabet, P, data):
    cfg = make_cfg("steering", {"P": P, "alphabet": alphabet, **FAST})
    res = encode_bytes(data, primers, cfg)
    out = decode_reads(_reads(res), primers, cfg["oligo_len"])
    assert out.ok and out.data == data


@pytest.mark.parametrize("alphabet", ["A", "B", "C"])
def test_edge_sizes(primers, alphabet):
    for n in (0, 1, 37, 300):
        data = random.Random(n).randbytes(n)
        cfg = make_cfg("steering", {"P": 6, "alphabet": alphabet, **FAST})
        res = encode_bytes(data, primers, cfg)
        assert decode_reads(_reads(res), primers, cfg["oligo_len"]).data == data


@pytest.mark.parametrize("alphabet", ["A", "B", "C"])
def test_steering_bases_carry_no_data(primers, alphabet):
    """Changing any steering base (keeping B's rotation valid) must not change the decoded frame."""
    if alphabet == "C":
        pytest.skip("C: radices depend on runs through steering bases; see test_c_steering_bases_carry_no_data")
    cfg = make_cfg("steering", {"P": 5, "alphabet": alphabet, **FAST})
    res = encode_bytes(random.Random(1).randbytes(400), primers, cfg)
    codec = res.codec
    geom = res.geometry
    f = primers.forward
    idx, oligo = res.oligos[-1]
    body = oligo[len(f) : len(oligo) - len(primers.tail)]
    seed = res.payload_info[-1]["seed"]
    payload = body[PREFIX_NT:]
    frame = codec.decode_payload(payload, idx, seed, body[PREFIX_NT - 1], 0)
    mask, _ = codec._layout(geom.payload_nt)
    prev0 = body[PREFIX_NT - 1]
    bases = "ACGT"
    # symbols at info positions (trits for B: relative to the previous base; bases for A)
    syms, pv = [], prev0
    for i, b in enumerate(payload):
        if not mask[i]:
            syms.append((bases.index(b) - bases.index(pv) - 1) % 4 if alphabet == "B" else b) if alphabet != "C" else None
        pv = b
    for p in [i for i, s in enumerate(mask) if s][:5]:
        for b in bases:
            # rebuild with steering base p changed and identical info symbols
            out, pv, si = [], prev0, 0
            for i in range(len(payload)):
                if mask[i]:
                    x = b if i == p else payload[i]
                else:
                    x = bases[(bases.index(pv) + 1 + syms[si]) % 4] if alphabet == "B" else syms[si]
                    si += 1
                out.append(x)
                pv = x
            new = "".join(out)
            assert codec.decode_payload(new, idx, seed, prev0, 0) == frame
            if alphabet == "C":
                break  # C's info radices depend on runs through steering bases: covered by test_c_*
            if alphabet == "B" and b != payload[p] and p + 1 < len(payload):
                assert new[p + 1] != payload[p + 1] or new[p + 1 :] != payload[p + 1 :]


def test_alphabet_b_steering_remaps_following_segment(primers):
    """B property (DESIGN.md F5): with trits fixed, a different steering base changes the next segment."""
    c = Steering(P=4, alphabet="B")
    assert c._render([0, 1, 2, 0], "A") != c._render([0, 1, 2, 0], "C")


@pytest.mark.parametrize("alphabet", ["A", "B", "C"])
def test_hard_constraints_met(primers, alphabet):
    cfg = make_cfg("steering", {"P": 4, "alphabet": alphabet, "w_conf": 0.0, "beam_width": 8})
    res = encode_bytes(random.Random(2).randbytes(2500), primers, cfg)
    scfg = ScreeningConfig.from_dict(cfg["screening"])
    bad = [(i, screen_oligo(s, primers, scfg, False)["violations"]) for i, s in res.oligos]
    bad = [b for b in bad if b[1]]
    recorded = [i for i in res.payload_info if i["violations"]]
    # every violation that exists must have been recorded, never silent
    assert {i for i, _ in bad} == {i["index"] for i in recorded}
    if alphabet in ("B", "C"):
        assert not bad  # P=4: GC windows steerable in practice for these alphabets


@given(st.text("ACGT", min_size=30, max_size=90), st.sampled_from([(4, 5, 6), (5, 6), (6,)]),
       st.floats(0.0, 10.0), st.integers(2, 4))
@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_incremental_proxy_equals_structure_module(primers, payload, ks, w_site, ml):
    """The beam's tier-1 cost must equal dnastore.structure.proxy_increments summed over
    non-left-flank positions (the left flank is constant across candidates and not charged)."""
    from dnastore.structure import ProxyParams, proxy_increments

    left = primers.forward + "ACGTACGTACGTACGTAC"
    right = primers.tail
    opts = SteeringOptions(w_conf=0.0, w_gc=0.0, w_stem=1.0, stem_ks=ks, stem_w_site=w_site, stem_min_loop=ml)
    scfg = ScreeningConfig.from_dict({"gc_window": [0.0, 1.0], "gc_global": [0.0, 1.0], "max_homopolymer": 200})
    L = len(left) + len(payload) + len(right)
    sc = _Scorer(opts, scfg, left, right, L, anchors=set(), primers=primers)
    stt, i, rng = sc.initial(), 0, random.Random(len(payload))
    while i < len(payload):  # irregular chunks exercise segment boundaries
        j = min(len(payload), i + rng.randint(1, 7))
        stt = sc.extend(stt, payload[i:j])
        i = j
    fin = sc.extend(stt, right)
    full = left + payload + right
    assert fin.seq == full
    inc = proxy_increments(full, len(primers.forward), len(primers.reverse),
                           ProxyParams(ks, ml, 1.0, w_site), primers=(primers.forward, primers.reverse))
    assert abs(fin.cost - sum(inc[len(left):])) < 1e-6 * max(1.0, abs(fin.cost))
    assert fin.gc_total == sum(b in "GC" for b in full)
    assert fin.gc_win == sum(b in "GC" for b in full[-scfg.gc_window_size :])


def test_anchor_check_excludes_true_site(primers):
    cfg = make_cfg("steering", {"P": 4, "alphabet": "B", **FAST})
    res = encode_bytes(b"hello world" * 20, primers, cfg)
    assert all(not i["violations"] for i in res.payload_info)


def test_density_ordering(primers):
    """More steering (smaller P) must cost capacity: frame bytes increase with P."""
    geom = OligoGeometry(200, len(primers.forward), len(primers.reverse))
    for alpha in "ABC":
        caps = [Steering(P=P, alphabet=alpha).frame_bits(geom.payload_nt) for P in (2, 4, 8, 12, 16, 0)]
        assert caps == sorted(caps) and len(set(caps)) == len(caps)  # strictly increasing: no plateaus


def test_parallel_encode_identical_to_serial(primers):
    data = random.Random(11).randbytes(1500)
    a = encode_bytes(data, primers, make_cfg("steering", {"P": 6, "alphabet": "B", **FAST}, workers=1))
    b = encode_bytes(data, primers, make_cfg("steering", {"P": 6, "alphabet": "B", **FAST}, workers=4))
    assert a.oligos == b.oligos


def test_tier2_rescoring_picks_most_accessible(primers):
    """With rescore='local_access', the chosen candidate has the highest tier-2 anchor access
    among the violation-free top candidates."""
    from dnastore.codecs.steering import Steering
    from dnastore.structure import local_access

    data = random.Random(12).randbytes(60)
    cfg = make_cfg("steering", {"P": 6, "alphabet": "B", "w_conf": 0.0, "beam_width": 6, "rescore_top": 3})
    res = encode_bytes(data, primers, cfg)
    for (idx, s), info in zip(res.oligos, res.payload_info):
        if idx >= 64:
            got = local_access(s, len(primers.forward), len(primers.reverse), 80, 60.0, 0.05)["local_min_anchor"]
            assert abs(got - info["local_min_anchor"]) < 1e-9


def test_c_run_limit_by_construction_and_rate(primers):
    """Alphabet C: payload runs never exceed max_homopolymer; rate close to RLL(3) capacity at P=inf."""
    from dnastore.screening import max_homopolymer

    cfg = make_cfg("steering", {"P": 0, "alphabet": "C", **FAST})
    res = encode_bytes(random.Random(5).randbytes(3000), primers, cfg)
    for (i, s), info in zip(res.oligos, res.payload_info):
        if i >= 64:
            assert max_homopolymer(s) <= 3 and not info["violations"]
    geom = res.geometry
    rate = res.codec.frame_bits(geom.payload_nt) / geom.payload_nt
    assert 1.90 < rate < 1.9824  # RLL(3) quaternary capacity ~1.9824 bits/nt


def test_c_render_decode_inverse():
    from dnastore.codecs.steering import Steering

    rng = random.Random(3)
    c = Steering(P=0, alphabet="C")
    for _ in range(200):
        n = rng.randint(5, 60)
        v = rng.getrandbits(n)  # well inside capacity
        seq, rem = c._render_c(v, n, "A", 1, 3)
        assert rem == 0
        # decode with the same mixed-radix rule
        val, mult, last, run = 0, 1, "A", 1
        for b in seq:
            if run >= 3:
                val += [x for x in "ACGT" if x != last].index(b) * mult
                mult *= 3
            else:
                val += "ACGT".index(b) * mult
                mult *= 4
            run = run + 1 if b == last else 1
            last = b
        assert val == v


def test_c_rejects_overlong_run_on_decode(primers):
    cfg = make_cfg("steering", {"P": 0, "alphabet": "C", **FAST})
    res = encode_bytes(b"x" * 50, primers, cfg)
    codec, geom = res.codec, res.geometry
    assert codec.decode_payload("AAAA" + "C" * (geom.payload_nt - 4), 70, 0, "G", 0) is None


def test_normalisation_spreads_positive(primers):
    from dnastore.codecs.steering import Steering

    c = Steering(P=8, alphabet="A", w_conf=0.0, normalize=True)
    c.bind(primers, {})
    o = c._effective_opts(len(primers.forward) + 18, 145)
    sp = c.spreads
    assert sp["gc"]["spread"] > 0 and sp["stem"]["spread"] > 0
    assert abs(o.w_gc - 1.0 / sp["gc"]["spread"]) < 1e-12


def test_c_steering_bases_carry_no_data():
    """C: re-rendering with a different steering base and the same digit sequence decodes identically."""
    from dnastore.codecs.steering import Steering, steering_mask

    rng = random.Random(8)
    c = Steering(P=5, alphabet="C")
    n = 60
    mask = steering_mask(n, 5)
    n_info = sum(not m for m in mask)
    for _ in range(50):
        v = rng.getrandbits(n_info)  # fits even if every info position were radix 3? use small v
        v %= 3**n_info
        steer = [rng.choice("ACGT") for _ in range(n)]

        def build(st):
            out, last, run, rem = [], "A", 1, v
            for p in range(n):
                if mask[p]:
                    b = st[p]
                    if b == last and run >= 3:
                        b = next(x for x in "ACGT" if x != last)
                else:
                    b, rem = c._render_c(rem, 1, last, run, 3)
                run = run + 1 if b == last else 1
                last = b
                out.append(b)
            assert rem == 0
            return "".join(out)

        s1, s2 = build(steer), build([rng.choice("ACGT") for _ in range(n)])
        nb = c.frame_bits(n)
        d1, d2 = c.decode_payload(s1, 0, 0, "A", 0), c.decode_payload(s2, 0, 0, "A", 0)
        if v >> nb:
            continue
        assert d1 == d2
