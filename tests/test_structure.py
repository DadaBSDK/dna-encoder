from __future__ import annotations

import random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dnastore.dnautil import revcomp
from dnastore.structure import (
    ProxyParams,
    StructureConfig,
    accessibility,
    mfe_at,
    proxy_increments,
    proxy_score,
    stem_weight,
    structure_report,
)

HAIRPIN = "GGGGCGCGCAAAAGCGCGCCCC"


def _random_oligo(primers, seed=0, n=160):
    rng = random.Random(seed)
    return primers.forward + "".join(rng.choice("ACGT") for _ in range(n)) + primers.tail


def test_temperature_and_salt_change_energies():
    e37 = mfe_at(HAIRPIN, 37.0, 0.05)
    e60 = mfe_at(HAIRPIN, 60.0, 0.05)
    e37_hi = mfe_at(HAIRPIN, 37.0, 1.021)
    assert e37 < e60 < 0  # hairpin is less stable at higher temperature
    assert e37_hi < e37  # higher salt stabilises


def test_accessibility_in_unit_interval(primers):
    acc = accessibility(_random_oligo(primers), len(primers.forward), len(primers.reverse), 60.0, 0.05)
    assert set(acc) >= {"site_access_R", "site_access_F", "anchor_access_R", "anchor_access_F",
                        "min_site_access", "min_anchor_access"}
    assert all(0.0 <= v <= 1.0 for v in acc.values())
    assert acc["min_anchor_access"] == min(acc["anchor_access_R"], acc["anchor_access_F"])


def test_sequestered_R_site_is_inaccessible(primers):
    base = _random_oligo(primers, seed=0)
    # Body ends with R itself + a 4-nt loop: R pairs with the adjacent rc(R) site -> hairpin.
    body = base[20:-20]
    trapped = primers.forward + body[: len(body) - 24] + primers.reverse + "TTTT" + primers.tail
    assert len(trapped) == len(base)
    a0 = accessibility(base, 20, 20, 60.0, 0.05)
    a1 = accessibility(trapped, 20, 20, 60.0, 0.05)
    assert a1["site_access_R"] < 0.2 < a0["site_access_R"]
    assert a1["anchor_access_R"] < 0.2
    assert proxy_score(trapped, 20, 20) > proxy_score(base, 20, 20)


def test_structure_report_keys_and_threshold(primers):
    cfg = StructureConfig(temperatures=(60.0,), min_anchor_access=0.0)
    rep = structure_report(_random_oligo(primers, 1), 20, 20, cfg)
    assert {"min_anchor_access_T60", "mfe_fwd_T60", "mfe_rc_T60"} <= set(rep)
    assert rep["access_ok"] is True


def test_structure_config_from_dict():
    c = StructureConfig.from_dict({"temperatures": [55, 37], "salt_molar": 0.1})
    assert c.temperatures == (55.0, 37.0) and c.salt_molar == 0.1
    with pytest.raises(ValueError):
        StructureConfig.from_dict({"bogus": 1})


def test_stem_weight_gc_heavier():
    assert stem_weight("GCGCG") > stem_weight("ATATA") > 0


@given(st.text(alphabet="ACGT", min_size=0, max_size=120), st.sampled_from([(4, 5, 6), (5, 6), (6, 7)]),
       st.integers(0, 5))
@settings(max_examples=60, deadline=None)
def test_proxy_sum_of_increments(seq, ks, loop):
    oligo = "GCTGATTTATGTGGCTTGCG" + seq + "TGTAGAGGGTTGAAGGCCAT"
    p = ProxyParams(ks=ks, min_loop=loop)
    inc = proxy_increments(oligo, 20, 20, p)
    assert len(inc) == len(oligo)
    assert proxy_score(oligo, 20, 20, p) == pytest.approx(sum(inc))
    # prefix property (what a left-to-right beam relies on): with the primers fixed,
    # increments of a prefix are unchanged by appending bases
    prim = ("GCTGATTTATGTGGCTTGCG", revcomp("TGTAGAGGGTTGAAGGCCAT"))
    assert proxy_increments(oligo, 20, 20, p, primers=prim) == inc
    longer = proxy_increments(oligo + "ACGT", 20, 20, p, primers=prim)
    assert all(a == pytest.approx(b) for a, b in zip(inc, longer))


def test_proxy_detects_planted_stem():
    flank = "GCTGATTTATGTGGCTTGCG"
    stem = "GACCTGAC"
    clean = flank + "A" * 40 + flank[::-1]
    planted = flank + stem + "TTTTT" + revcomp(stem) + "A" * 19 + flank[::-1]
    p = ProxyParams(ks=(5, 6), w_site=0.0)
    assert proxy_score(planted, 20, 20, p) > proxy_score(clean, 20, 20, p)
