"""Noisy end-to-end: channel -> decoder recovers the file; decoder reports are honest."""

from __future__ import annotations

import random

import numpy as np
import pytest

from dnastore.channel import ChannelParams, simulate
from dnastore.decoder import decode_reads
from dnastore.encoder import encode_bytes

from .conftest import make_cfg

FAST = {"w_conf": 0.0, "beam_width": 4, "rescore": None, "normalize": False}


@pytest.mark.parametrize("name,params", [("naive2bit", {"whiten": True}), ("goldman", {}),
                                         ("steering", {"P": 8, "alphabet": "A", **FAST}),
                                         ("steering", {"P": 8, "alphabet": "C", **FAST})])
def test_noisy_channel_recovery(primers, name, params):
    data = random.Random(1).randbytes(3000)
    cfg = make_cfg(name, params)
    res = encode_bytes(data, primers, cfg)
    oligos = [s for _, s in res.oligos]
    ch = ChannelParams.from_dict({"p_sub": 0.01, "p_ins": 0.005, "p_del": 0.005, "dropout": 0.02,
                                  "mean_coverage": 12, "coverage": "lognormal", "lognormal_sigma": 0.3})
    reads, truth = simulate(oligos, ch, np.random.default_rng(0))
    out = decode_reads([r.seq for r in reads], primers, cfg["oligo_len"], workers=4)
    assert out.ok and out.data == data, out.report
    assert out.report["reverse_strand"] > 0.3 * out.report["reads"]


def test_insufficient_coverage_fails_loudly(primers):
    data = random.Random(2).randbytes(3000)
    cfg = make_cfg("goldman", {})
    res = encode_bytes(data, primers, cfg)
    ch = ChannelParams.from_dict({"p_sub": 0.03, "p_ins": 0.03, "p_del": 0.03, "mean_coverage": 1.5})
    reads, _ = simulate([s for _, s in res.oligos], ch, np.random.default_rng(1))
    out = decode_reads([r.seq for r in reads], primers, cfg["oligo_len"])
    assert not out.ok and out.data is None and "error" in out.report
