"""End-to-end smoke test of the squigulator -> Dorado channel (skipped if tools are missing)."""

from __future__ import annotations

import numpy as np
import pytest

from dnastore import nanopore
from dnastore.dnautil import revcomp

pytestmark = [pytest.mark.slow, pytest.mark.skipif(not nanopore.tools_available("hac"), reason="nanopore tools not installed")]


def test_write_reference_strand_and_ids(tmp_path):
    truth = nanopore.write_reference(["ACGTACGTAA", "GGGCCCAATT"], [3, 2], tmp_path / "r.fa", seed=1, flanks=False)
    seqs = [l.strip() for l in open(tmp_path / "r.fa") if not l.startswith(">")]
    assert len(truth) == 5 and truth[4]["read_id"].endswith("000000000005")
    for t, s in zip(truth, seqs):
        src = ["ACGTACGTAA", "GGGCCCAATT"][t["oligo"]]
        assert s == (revcomp(src) if t["reverse"] else src)


def test_end_to_end_small(tmp_path, primers):
    import edlib

    rng = np.random.default_rng(0)
    oligos = [primers.forward + "".join(rng.choice(list("ACGT"), 160)) + primers.tail for _ in range(10)]
    fq, meta = nanopore.run_pipeline(oligos, [3] * 10, tmp_path, seed=7, model="hac")
    reads = nanopore.read_fastq(fq)
    assert reads
    truth = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(tmp_path / "truth.tsv").readlines()[1:]}
    best = 0.0
    for rid, seq, _ in reads:
        src = oligos[truth[rid]]
        for s in (src, revcomp(src)):
            r = edlib.align(s, seq, mode="HW", task="distance")
            best = max(best, 1 - r["editDistance"] / len(s))
    assert best >= 0.8
