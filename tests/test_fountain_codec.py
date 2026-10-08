"""Fountain arm: screening, round trips, and loud failure."""

from __future__ import annotations

import random

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from dnastore.decoder import decode_reads
from dnastore.dnautil import revcomp
from dnastore.encoder import encode_bytes
from dnastore.screening import ScreeningConfig, screen_oligo

from .conftest import make_cfg


def _reads(res, drop=frozenset(), seed=0):
    rng = random.Random(seed)
    reads = [s if rng.random() < 0.5 else revcomp(s) for i, s in res.oligos if i not in drop]
    rng.shuffle(reads)
    return reads


@given(data=st.binary(max_size=2500))
@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_roundtrip_noiseless(primers, data):
    cfg = make_cfg("fountain", {"overhead": 0.6})
    res = encode_bytes(data, primers, cfg)
    out = decode_reads(_reads(res), primers, cfg["oligo_len"])
    if res.codec.k < 10:
        # LT with a handful of segments fails with non-negligible probability even noiselessly
        # (e.g. k = 2, 4 droplets); it must then fail loudly, never return wrong data.
        assert (out.ok and out.data == data) or (not out.ok and out.data is None), out.report
    else:
        assert out.ok and out.data == data, out.report


def test_every_emitted_droplet_passes_screening(primers):
    cfg = make_cfg("fountain", {"overhead": 0.2})
    res = encode_bytes(random.Random(1).randbytes(4000), primers, cfg)
    sc = ScreeningConfig.from_dict(cfg["screening"])
    assert all(not screen_oligo(s, primers, sc, False)["violations"] for i, s in res.oligos if i >= 64)
    assert res.layout.rejected > 0  # whitened 2-bit droplets fail screening often (DESIGN.md F4)


def test_survives_droplet_loss(primers):
    data = random.Random(2).randbytes(8000)
    cfg = make_cfg("fountain", {"overhead": 0.25})
    res = encode_bytes(data, primers, cfg)
    data_idx = [i for i, _ in res.oligos if i >= 64]
    drop = set(random.Random(3).sample(data_idx, len(data_idx) // 10))
    out = decode_reads(_reads(res, drop), primers, cfg["oligo_len"])
    assert out.ok and out.data == data


def test_too_few_droplets_fails_loudly(primers):
    data = random.Random(4).randbytes(4000)
    cfg = make_cfg("fountain", {"overhead": 0.05})
    res = encode_bytes(data, primers, cfg)
    data_idx = [i for i, _ in res.oligos if i >= 64]
    drop = set(data_idx[: len(data_idx) // 5])
    out = decode_reads(_reads(res, drop), primers, cfg["oligo_len"])
    assert not out.ok and out.data is None and out.report["error"] == "LT decode failure"
