"""Noiseless round-trip tests: every codec must return byte-identical files."""

from __future__ import annotations

import hashlib
import json
import random

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from typer.testing import CliRunner

from dnastore.cli import app
from dnastore.decoder import decode_reads
from dnastore.dnautil import revcomp
from dnastore.ecc import CRC_BITS
from dnastore.encoder import encode_bytes
from dnastore.oligo import HEADER_INDEX_SPACE, INDEX_TRITS, PREFIX_NT, decode_index
from dnastore.screening import max_homopolymer

from .conftest import CODEC_CONFIGS, ROOT, make_cfg


def _pool_reads(res, seed=0, drop=()):
    rng = random.Random(seed)
    reads = [s if rng.random() < 0.5 else revcomp(s) for i, s in res.oligos if i not in drop]
    rng.shuffle(reads)
    return reads


def _roundtrip(data, name, params, primers, **kw):
    cfg = make_cfg(name, params, **kw)
    res = encode_bytes(data, primers, cfg)
    assert all(len(s) == cfg["oligo_len"] for _, s in res.oligos)
    out = decode_reads(_pool_reads(res), primers, cfg["oligo_len"])
    assert out.ok, out.report
    assert hashlib.sha256(out.data).digest() == hashlib.sha256(data).digest()
    return res, out


@pytest.mark.parametrize("name,params", CODEC_CONFIGS)
@given(data=st.binary(max_size=3000))
@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_random_bytes_roundtrip(primers, name, params, data):
    _roundtrip(data, name, params, primers)


@pytest.mark.parametrize("name,params", CODEC_CONFIGS)
def test_edge_sizes(primers, name, params):
    cfg = make_cfg(name, params)
    from dnastore.codecs import make_codec
    from dnastore.ecc import default_k_max
    from dnastore.oligo import OligoGeometry

    geom = OligoGeometry(cfg["oligo_len"], len(primers.forward), len(primers.reverse))
    b = -(-(make_codec(name, params).frame_bits(geom.payload_nt) - CRC_BITS) // 8)  # RS row bytes
    kmax = default_k_max(cfg["ecc"]["parity_permille"])
    rng = random.Random(7)
    for n in [0, 1, b - 1, b, b + 1, kmax * b - 1, kmax * b, kmax * b + 1]:
        data = bytes(rng.getrandbits(8) for _ in range(n))
        res, out = _roundtrip(data, name, params, primers)
        if n == kmax * b + 1:
            assert len(res.layout.blocks) == 2


@pytest.mark.parametrize("name,params", CODEC_CONFIGS)
def test_compressible_text_and_extension(primers, name, params):
    text = (ROOT / "docs" / "DESIGN.md").read_bytes()
    res, out = _roundtrip(text, name, params, primers)
    assert res.info.compressed and out.extension == ".txt"


def test_survives_dropout_within_parity_and_header_copies(primers):
    data = random.Random(3).randbytes(4000)
    cfg = make_cfg("goldman", {})
    res = encode_bytes(data, primers, cfg)
    start, (k, p) = res.layout.block_starts()[0], res.layout.blocks[0]
    drop = {HEADER_INDEX_SPACE + i for i in range(p - 2)}  # bit-packed: allow for row-straddling bytes
    drop |= {i for i, _ in res.oligos if i < HEADER_INDEX_SPACE and i % cfg["header_copies"] != 0}  # all but 1 header copy
    out = decode_reads(_pool_reads(res, drop=drop), primers, cfg["oligo_len"])
    assert out.ok and out.data == data


def test_fails_loudly_beyond_parity(primers):
    data = random.Random(4).randbytes(4000)
    cfg = make_cfg("goldman", {})
    res = encode_bytes(data, primers, cfg)
    k, p = res.layout.blocks[0]
    drop = {HEADER_INDEX_SPACE + i for i in range(p + 2)}
    out = decode_reads(_pool_reads(res, drop=drop), primers, cfg["oligo_len"])
    assert not out.ok and out.data is None and out.report["error"] == "RS failure"


def test_missing_header_fails_loudly(primers):
    cfg = make_cfg("goldman", {})
    res = encode_bytes(b"abc" * 100, primers, cfg)
    out = decode_reads(_pool_reads(res, drop=set(range(HEADER_INDEX_SPACE))), primers, cfg["oligo_len"])
    assert not out.ok and "header" in out.report["error"]


def test_header_copies_are_distinct_sequences(primers):
    res = encode_bytes(b"hello", primers, make_cfg("goldman", {}))
    hdr = [s for i, s in res.oligos if i < HEADER_INDEX_SPACE]
    assert len(set(hdr)) == len(hdr) == res.n_header_oligos


def test_index_field_and_ternary_payload_homopolymer_free(primers):
    res = encode_bytes(random.Random(5).randbytes(2000), primers, make_cfg("goldman", {}))
    f = primers.forward
    for idx, s in res.oligos:
        body = s[len(f) : len(s) - len(primers.tail)]
        assert decode_index(body[:INDEX_TRITS], f[-1]) == idx
        assert max_homopolymer(f[-1] + body) == 1  # rotating code end to end (index, seed, payload)


def test_cli_roundtrip(tmp_path):
    src = tmp_path / "sample.txt"
    src.write_bytes(b"DNA storage CLI round trip\n" * 50)
    runner = CliRunner()
    r = runner.invoke(app, ["encode", str(src), "--out-dir", str(tmp_path), "--config", str(ROOT / "configs" / "default.yaml")])
    assert r.exit_code == 0, r.output
    manifest = json.loads((tmp_path / "sample.txt.goldman.manifest.json").read_text())
    assert manifest["metrics"]["n_oligos"] > 0
    r = runner.invoke(app, ["decode", str(tmp_path / "sample.txt.goldman.fasta"), "--out", str(tmp_path / "rec"),
                            "--config", str(ROOT / "configs" / "default.yaml")])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "rec.txt").read_bytes() == src.read_bytes()


@pytest.mark.parametrize("name,params", CODEC_CONFIGS)
def test_primer_anchor_is_hard_for_all_arms(primers, name, params):
    """3'-anchored primer matches are rerolled away in every arm (DESIGN.md Part 0, rev. 2)."""
    from dnastore.encoder import framework_violations
    from dnastore.oligo import decode_prefix

    cfg = make_cfg(name, params)
    # Bytes whose raw 2-bit mapping spells F's 3'-terminal 8-mer, forcing hits in naive2bit.
    from dnastore.dnautil import bases_to_bytes_2bit

    data = bases_to_bytes_2bit(primers.forward[-8:]) * 3000
    res = encode_bytes(data + random.Random(9).randbytes(6000), primers, cfg)
    reroll = 0
    for (idx, s), info in zip(res.oligos, res.payload_info):
        assert info["index"] == idx
        assert framework_violations(s, primers, cfg["screening"]) == [] or info["violations"] == ["primer_anchor"]
        assert decode_prefix(s[len(primers.forward):], primers.forward[-1])[1] == info["seed"]
        reroll += info["seed"] > 0
    assert all(not i["violations"] for i in res.payload_info)
    if name == "naive2bit" and not params:
        assert reroll > 0  # the planted motif must have forced rerolls
    out = decode_reads(_pool_reads(res), primers, cfg["oligo_len"])
    assert out.ok
