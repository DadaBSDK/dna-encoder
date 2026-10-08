import random

import pytest

from dnastore.dnautil import revcomp
from dnastore.primers import PrimerPair
from dnastore.screening import (
    ScreeningConfig,
    gc_fraction,
    gc_window_extremes,
    max_homopolymer,
    mfe,
    mfe_both,
    primer_hits,
    screen_oligo,
    screen_pool,
)

PRIMERS = PrimerPair("GCTGATTTATGTGGCTTGCG", "ATGGCCTTCAACCCTCTACA")


def rand_seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def body_without_hits(n, seed):
    """Random body with no spurious primer hits when wrapped by PRIMERS."""
    for s in range(seed, seed + 1000):
        b = rand_seq(n, s)
        if primer_hits(PRIMERS.forward + b + PRIMERS.tail, PRIMERS, 3, 8, 0)["full"] == 0:
            ph = primer_hits(PRIMERS.forward + b + PRIMERS.tail, PRIMERS, 3, 8, 0)
            if ph["anchor"] == 0:
                return b
    raise AssertionError


def test_gc_and_homopolymer():
    assert gc_fraction("") == 0.0
    assert gc_fraction("GCAT") == 0.5
    assert max_homopolymer("") == 0
    assert max_homopolymer("A") == 1
    assert max_homopolymer("ACGTTTTAC") == 4
    assert max_homopolymer("AAACCCGGGTTT") == 3


def test_gc_window():
    seq = "A" * 10 + "G" * 10 + "A" * 10
    lo, hi = gc_window_extremes(seq, 10)
    assert lo == 0.0 and hi == 1.0
    assert gc_window_extremes("GCAT", 20) == (0.5, 0.5)
    lo, hi = gc_window_extremes("GCGCATAT", 4)
    assert (lo, hi) == (0.0, 1.0)


def test_true_sites_not_counted():
    oligo = PRIMERS.forward + body_without_hits(150, 0) + PRIMERS.tail
    ph = primer_hits(oligo, PRIMERS, 3, 8, 0)
    assert ph["full"] == 0 and ph["anchor"] == 0


def _mutate(s, k, seed):
    r = random.Random(seed)
    s = list(s)
    for pos in r.sample(range(len(s)), k):
        s[pos] = r.choice([b for b in "ACGT" if b != s[pos]])
    return "".join(s)


@pytest.mark.parametrize("mm", [0, 1, 2, 3])
@pytest.mark.parametrize("which", ["F", "R", "rcF", "rcR"])
def test_injected_primer_copies_detected(mm, which):
    body = body_without_hits(150, 10)
    q = PRIMERS.orientations()[which]
    # mutate only the 5' half of the primer orientation string so anchored check is independent
    ins = _mutate(q, mm, seed=mm)
    oligo = PRIMERS.forward + body[:60] + ins + body[60 + len(ins):] + PRIMERS.tail
    ph = primer_hits(oligo, PRIMERS, 3, 8, 0)
    assert any(name == which and pos == len(PRIMERS.forward) + 60 for name, pos, _ in ph["full_positions"])
    # 4 mismatches are beyond the threshold
    ins4 = _mutate(q, 4, seed=99)
    oligo4 = PRIMERS.forward + body[:60] + ins4 + body[60 + len(ins4):] + PRIMERS.tail
    ph4 = primer_hits(oligo4, PRIMERS, 3, 8, 0)
    assert not any(name == which and pos == len(PRIMERS.forward) + 60 for name, pos, _ in ph4["full_positions"])


def test_anchor_hits_both_orientations():
    body = body_without_hits(150, 20)
    f3 = PRIMERS.forward[-8:]
    rc_f3 = revcomp(f3)
    for frag, name in ((f3, "F3"), (rc_f3, "rcF3")):
        oligo = PRIMERS.forward + body[:50] + frag + body[58:] + PRIMERS.tail
        ph = primer_hits(oligo, PRIMERS, 3, 8, 0)
        assert (name, len(PRIMERS.forward) + 50, 0) in ph["anchor_positions"]


def test_mfe_hairpin_vs_random():
    stem = "GCGCGGCCGCATGC"
    hairpin = stem + "TTTT" + revcomp(stem)
    assert mfe(hairpin) < -15
    rnd = rand_seq(len(hairpin), 3)
    assert mfe(hairpin) < mfe(rnd)
    f, r = mfe_both(hairpin)
    assert f < -15 and r < -15


def test_screen_oligo_violations():
    cfg = ScreeningConfig()
    good = PRIMERS.forward + body_without_hits(150, 30) + PRIMERS.tail
    res = screen_oligo(good, PRIMERS, cfg, compute_mfe=False)
    assert res["mfe_fwd"] is None
    bad = PRIMERS.forward + "A" * 30 + body_without_hits(120, 40) + PRIMERS.tail
    res = screen_oligo(bad, PRIMERS, cfg, compute_mfe=False)
    assert "homopolymer" in res["violations"]
    assert "gc_window" in res["violations"]
    # mfe constraint only active with threshold and in hard list
    cfg2 = ScreeningConfig(mfe_threshold=0.0, hard=("mfe",))
    res = screen_oligo(good, PRIMERS, cfg2, compute_mfe=True)
    assert res["violations"] == ["mfe"]  # any structured oligo has MFE < 0
    assert res["mfe_fwd"] is not None and res["mfe_rc"] is not None


def test_config_from_dict():
    c = ScreeningConfig.from_dict({"gc_global": [0.3, 0.7], "hard": ["homopolymer"]})
    assert c.gc_global == (0.3, 0.7) and c.hard == ("homopolymer",)
    with pytest.raises(ValueError):
        ScreeningConfig.from_dict({"nope": 1})


def test_screen_pool_parallel_matches_serial():
    cfg = ScreeningConfig()
    pool = [PRIMERS.forward + rand_seq(60, s) + PRIMERS.tail for s in range(10)]
    a = screen_pool(pool, PRIMERS, cfg, compute_mfe=True, workers=2)
    b = screen_pool(pool, PRIMERS, cfg, compute_mfe=True, workers=1)
    assert a == b
