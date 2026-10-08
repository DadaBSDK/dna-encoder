from __future__ import annotations

import numpy as np

from dnastore.channel import ChannelParams, simulate
from dnastore.consensus import cluster, consensus, locate_and_orient

F = "GCTGATTTATGTGGCTTGCG"
TAIL = "TGTAGAGGGTTGAAGGCCAT"  # revcomp of R = ATGGCCTTCAACCCTCTACA


def _rand(n, rng):
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_locate_and_orient_both_strands_with_primer_errors():
    rng = np.random.default_rng(0)
    body = _rand(160, rng)
    oligo = F + body + TAIL
    got = locate_and_orient(oligo, F, TAIL)
    assert got[0] == body and got[1] is False
    got = locate_and_orient(_rc(oligo), F, TAIL)
    assert got[0] == body and got[1] is True
    # 2 substitutions + 1 deletion in F, 1 insertion in tail, reverse strand
    f_bad = "GCAGATTTATGTGCTTGCG"
    f_bad = f_bad[:5] + "C" + f_bad[6:]
    t_bad = TAIL[:7] + "A" + TAIL[7:]
    got = locate_and_orient(_rc(f_bad + body + t_bad), F, TAIL)
    assert got is not None and got[1] is True and got[0] == body
    assert locate_and_orient(_rand(200, rng), F, TAIL) is None


def test_locate_and_orient_with_adapter_flanks():
    # untrimmed nanopore reads: leader/adapter bases before F and after rc(R), either strand
    from dnastore.nanopore import FLANK_3, FLANK_5

    rng = np.random.default_rng(3)
    body = _rand(160, rng)
    oligo = F + body + TAIL
    got = locate_and_orient(FLANK_5 + oligo + FLANK_3, F, TAIL)
    assert got is not None and got[0] == body and got[1] is False
    got = locate_and_orient(FLANK_5 + _rc(oligo) + FLANK_3, F, TAIL)
    assert got is not None and got[0] == body and got[1] is True
    # no false hits on random reads of flanked length
    assert all(locate_and_orient(_rand(264, rng), F, TAIL) is None for _ in range(200))


def _accuracy(rate, cov, trials=60, seed=1):
    rng = np.random.default_rng(seed)
    ok = 0
    for t in range(trials):
        body = _rand(160, rng)
        p = ChannelParams(p_sub=rate / 3, p_ins=rate / 3, p_del=rate / 3, mean_coverage=cov, rc_frac=0.0, seed=seed * 1000 + t)
        reads, _ = simulate([body], p)
        seq, _ = consensus([r.seq for r in reads], 160, rng)
        ok += seq == body
    return ok / trials


def test_consensus_recovers_body_at_3pct_ids_10x():
    assert _accuracy(0.03, 10) >= 0.9


def test_consensus_monotone_in_coverage_at_3pct():
    """Regression: more reads must not hurt (column-majority-only consensus did)."""
    assert _accuracy(0.03, 20, trials=60, seed=3) >= 0.95


def test_consensus_exact_path_and_single():
    rng = np.random.default_rng(2)
    body = _rand(100, rng)
    reads = [body] * 3 + [body[:10] + "A" + body[11:]] * 2
    seq, info = consensus(reads + [body] * 2, 100, rng)
    assert seq == body and info["method"] == "align_vote(colmaj_init)"
    seq, info = consensus([body], 100)
    assert seq == body and info["method"] == "single"


def test_cluster_separates_distinct_bodies():
    rng = np.random.default_rng(3)
    a, b = _rand(150, rng), _rand(150, rng)
    p = ChannelParams(p_sub=0.01, p_ins=0.01, p_del=0.01, mean_coverage=8, rc_frac=0.0, seed=4)
    ra, _ = simulate([a], p)
    rb, _ = simulate([b], ChannelParams(**{**p.__dict__, "seed": 5}))
    bodies = [r.seq for r in ra] + [r.seq for r in rb]
    cl = cluster(bodies, 0.15, rng)
    assert len(cl) == 2
    assert {frozenset(c) for c in cl} == {frozenset(range(len(ra))), frozenset(range(len(ra), len(bodies)))}
