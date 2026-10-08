"""CRC-guided homopolymer repair and the two decoder strengths."""

from __future__ import annotations

import random

import numpy as np

from dnastore.decoder import decode_reads
from dnastore.dnautil import revcomp
from dnastore.encoder import encode_bytes
from dnastore.oligo import HEADER_INDEX_SPACE
from dnastore.repair import candidates, repair, run_support, runs

from .conftest import make_cfg


def test_candidates_are_runlength_edits_of_target_length():
    ref = "AACGTTTGCA"
    for target in (9, 10, 11, 12, 8):
        out = [c for c, _ in candidates(ref, target, budget=10**6)]
        assert out and all(len(c) == target for c in out)
        assert len(set(out)) == len(out)
        for c in out:  # every candidate keeps the run order of the reference (run bases)
            assert [b for b, _ in runs(c)] == [b for b, _ in runs(ref)] or len(runs(c)) < len(runs(ref))
    assert list(candidates(ref, 13)) == []  # more than 2 edits away


def test_support_ranks_the_true_edit_first():
    rng = np.random.default_rng(0)
    true = "".join(rng.choice(list("ACGT"), 120))
    i = next(j for j in range(40, 120) if true[j] == true[j + 1])  # a run of length >= 2
    bad = true[:i] + true[i + 1 :]  # systematic deletion inside that run
    reads = [bad] * 6 + [true] * 4
    sup = run_support(reads, bad)
    got = repair(bad, reads, len(true), lambda c: c == true)
    assert got.body == true and got.tried <= 3
    assert max(p for _, p in sup) >= 0.4


def _reads(res, mutate: dict[int, str], copies=6, seed=0):
    rng = random.Random(seed)
    out = []
    for idx, s in res.oligos:
        s = mutate.get(idx, s)
        for _ in range(copies):
            out.append(s if rng.random() < 0.5 else revcomp(s))
    rng.shuffle(out)
    return out


def test_repair_strength_recovers_systematic_homopolymer_errors(primers):
    cfg = make_cfg("naive2bit", {"whiten": True}, ecc={"parity_permille": 50, "k_max": None})
    res = encode_bytes(random.Random(5).randbytes(3000), primers, cfg)
    lf, lt = len(primers.forward), len(primers.tail)
    data_oligos = [(i, s) for i, s in res.oligos if i >= HEADER_INDEX_SPACE]
    mutate = {}
    for idx, s in data_oligos[: len(data_oligos) // 3]:  # far more than the 5% parity can erase
        body = s[lf:-lt]
        rr = runs(body)
        k = max(range(5, len(rr)), key=lambda j: rr[j][1])  # longest run after the index
        e = [n for _, n in rr]
        e[k] -= 1  # every read of this oligo carries the same run-length error
        body2 = "".join(b * n for (b, _), n in zip(rr, e))
        mutate[idx] = s[:lf] + body2 + s[-lt:]
    reads = _reads(res, mutate)
    weak = decode_reads(reads, primers, cfg["oligo_len"], strength="erasure")
    strong = decode_reads(reads, primers, cfg["oligo_len"], strength="repair")
    assert not weak.ok
    assert strong.ok and strong.data == random.Random(5).randbytes(3000), strong.report
    assert strong.report["repaired"] == len(mutate)
