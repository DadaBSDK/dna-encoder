from pathlib import Path

from dnastore.dnautil import revcomp
from dnastore.primers import PrimerPair, design_primer_pair, load_primers, save_primers, validate_primer_pair

ROOT = Path(__file__).resolve().parents[1]


def test_design_deterministic_and_valid():
    p1, r1 = design_primer_pair(1)
    p2, _ = design_primer_pair(1)
    assert p1 == p2
    assert r1["ok"], r1["failures"]
    assert len(p1.forward) == len(p1.reverse) == 20


def test_shipped_primers_validate():
    pair = load_primers(ROOT / "configs" / "primers.yaml")
    rep = validate_primer_pair(pair)
    assert rep["ok"], rep["failures"]


def test_tail_is_revcomp():
    p = PrimerPair("ACGTACGTAC", "GGGTTTCCCA")
    assert p.tail == revcomp("GGGTTTCCCA")


def test_validation_rejects_bad_pairs():
    # Homopolymer, extreme GC, palindromic.
    bad = PrimerPair("AAAAAAAAAACCCCCCCCCC", "GCGCGCGCGCGCGCGCGCGC")
    rep = validate_primer_pair(bad)
    assert not rep["ok"]
    assert any("homopolymer" in f for f in rep["failures"])
    # Identical primers fail the cross-similarity check.
    good, _ = design_primer_pair(1)
    same = PrimerPair(good.forward, good.forward)
    assert any("cross" in f for f in validate_primer_pair(same)["failures"])


def test_save_load_roundtrip(tmp_path):
    p, r = design_primer_pair(3)
    save_primers(p, tmp_path / "p.yaml", r)
    assert load_primers(tmp_path / "p.yaml") == p
