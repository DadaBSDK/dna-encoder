"""Adversarial tests: hostile inputs, corruption, contamination and misuse.

Two kinds of test live here:

* **Guards** pin down safety properties that hold today (never return wrong bytes, reject
  corrupt archives cleanly, RS recovers exactly up to capacity, web isolation).
* **Known bugs** are ``xfail(strict=True)`` with a ``BUG:`` reason. They document a
  reproducible defect without turning CI red; once the defect is fixed the test XPASSes,
  strict mode fails the run, and the marker should be removed.
"""

from __future__ import annotations

import gzip
import http.client
import json
import random
import threading
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pytest
from typer.testing import CliRunner

from dnastore.cli import app
from dnastore.config import load_config
from dnastore.decoder import decode_reads, read_fasta
from dnastore.dnautil import revcomp
from dnastore.ecc import _ee_column, parity_rows, plan_layout, rs_decode, rs_encode, rs_encode_matrix
from dnastore.encoder import encode_bytes
from dnastore.stream import decode_file, encode_file

from .conftest import ROOT, make_cfg


def _mutate(s: str, rng: random.Random, rate: float) -> str:
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3:
            continue
        if r < 2 * rate / 3:
            out.append(rng.choice("ACGT"))
        elif r < rate:
            c = rng.choice("ACGT".replace(c, ""))
        out.append(c)
    return "".join(out)


# --------------------------------------------------------------------------- guards: decoder safety


@pytest.mark.parametrize("name,params", [("goldman", {}), ("naive2bit", {"whiten": True}),
                                         ("fountain", {"overhead": 0.3})])
def test_corrupted_reads_never_crash_or_return_wrong_bytes(primers, name, params):
    """Heavy IDS noise, truncation, chimeras, junk and dropout: fail cleanly or succeed exactly."""
    rng = random.Random(hash(name) & 0xFFFF)
    data = rng.randbytes(700)
    cfg = make_cfg(name, params, global_seed=rng.randrange(2**32))
    res = encode_bytes(data, primers, cfg)
    for trial in range(6):
        rate, depth, drop = rng.choice([0.0, 0.03, 0.15]), rng.choice([1, 4]), rng.choice([0.0, 0.3])
        reads = []
        for _, s in res.oligos:
            if rng.random() < drop:
                continue
            for _ in range(rng.randint(1, depth)):
                r = _mutate(s, rng, rate)
                r = revcomp(r) if rng.random() < 0.5 else r
                if rng.random() < 0.05:
                    r = r[: rng.randrange(len(r) + 1)]
                if rng.random() < 0.03:
                    r = "".join(rng.choice("ACGT") for _ in range(rng.randrange(400)))
                if rng.random() < 0.03:
                    r += r
                reads.append(r)
        rng.shuffle(reads)
        out = decode_reads(reads, primers, cfg["oligo_len"], strength=rng.choice(["erasure", "repair"]),
                           grouping=rng.choice(["cluster", "index"]))
        assert not out.ok or out.data == data, (trial, rate, depth, drop)


def test_mixed_pools_never_return_wrong_file(primers):
    """Two files encoded with the same primers and indices, mixed at several ratios."""
    a, b = random.Random(1).randbytes(1500), random.Random(2).randbytes(1500)
    ra = [s for _, s in encode_bytes(a, primers, make_cfg("goldman", {})).oligos]
    rb = [s for _, s in encode_bytes(b, primers, make_cfg("goldman", {})).oligos]
    for reads in (ra * 3 + rb * 3, ra * 5 + rb, ra + rb * 5):
        for grouping in ("cluster", "index"):
            out = decode_reads(reads, primers, 200, grouping=grouping)
            assert not out.ok or out.data in (a, b)


def test_decoder_ignores_garbage_only_input(primers):
    rng = random.Random(0)
    junk = ["".join(rng.choice("ACGTN") for _ in range(rng.randrange(0, 300))) for _ in range(200)]
    out = decode_reads(junk + ["", "N" * 200], primers, 200)
    assert not out.ok and out.report["error"] == "header incomplete"


# --------------------------------------------------------------------------- guards: outer RS


def test_rs_recovers_any_erasure_pattern_within_capacity():
    """Drop random oligos (bit-packed, so they straddle rows); whenever every block column has
    at most p erasures the data must come back, and it must never come back wrong."""
    rng = random.Random(0)
    for _ in range(150):
        pb = rng.choice([8, 9, 17, 64, 65, 333])
        perm = rng.choice([0, 1, 150, 1000, 3000])
        kmax = rng.choice([None, 1, 7, 50])
        if kmax and kmax + parity_rows(kmax, perm) > 255:
            continue
        stored = rng.randbytes(rng.choice([0, 1, 100, 2000]))
        lay = plan_layout(len(stored), pb, perm, kmax)
        frames = rs_encode(stored, lay)
        drop = {q for q in range(lay.n_oligos) if rng.random() < 0.15}
        have = np.ones(lay.n_oligos * pb, bool)
        for q in drop:
            have[q * pb:(q + 1) * pb] = False
        within = True
        if lay.n_rows:
            known = have[:lay.stream_bits].reshape(-1, 8).all(1).reshape(lay.n_rows, lay.row_bytes)
            within = all((~known[s:s + k + p]).sum(0).max() <= p for s, (k, p) in zip(lay.block_starts(), lay.blocks))
        got = rs_decode({q: f for q, f in enumerate(frames) if q not in drop}, lay)
        assert got.data is None or got.data == stored
        if within:
            assert got.data == stored


def test_rs_single_undetected_bad_oligo_never_silent():
    rng = random.Random(3)
    lay = plan_layout(3000, 64, 150)
    for _ in range(40):
        stored = rng.randbytes(3000)
        rec = dict(enumerate(rs_encode(stored, lay)))
        q = rng.randrange(lay.n_oligos)
        rec[q] ^= 1 << rng.randrange(64)
        for d in rng.sample([i for i in rec if i != q], rng.randint(0, 3)):
            del rec[d]
        got = rs_decode(rec, lay)
        assert got.data is None or got.data == stored


# --------------------------------------------------------------------------- guards: streaming archives


def test_stream_archive_corruption_always_rejected(tmp_path):
    rng = random.Random(0)
    src = tmp_path / "x.bin"
    src.write_bytes(rng.randbytes(2000))
    encode_file(src, tmp_path / "x.dna")
    blob = (tmp_path / "x.dna").read_bytes()
    variants = [blob[:rng.randrange(len(blob))] for _ in range(40)]
    for _ in range(80):
        b = bytearray(blob)
        b[rng.randrange(len(b))] ^= 1 << rng.randrange(8)
        variants.append(bytes(b))
    variants += [blob.replace(b"\n", b"\r\n"), blob.lower()]
    for n, v in enumerate(variants):
        if v == blob:
            continue
        (tmp_path / "m.dna").write_bytes(v)
        out = tmp_path / f"out{n}"
        with pytest.raises(ValueError):
            decode_file(tmp_path / "m.dna", out)
        assert not out.exists()
    assert not [p for p in tmp_path.iterdir() if p.name.startswith(".dnastore-")]


# --------------------------------------------------------------------------- guards: web isolation


@pytest.fixture
def workbench(tmp_path):
    from dnastore.web import make_server

    srv = make_server(port=0, root=tmp_path)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv
    srv.shutdown()
    srv.server_close()
    srv.workbench.executor.shutdown(wait=True)


def _req(port, method, path, body=None, headers=None, host=None):
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
    c.putrequest(method, path, skip_host=True, skip_accept_encoding=True)
    c.putheader("Host", host or f"127.0.0.1:{port}")
    for k, v in (headers or {}).items():
        c.putheader(k, v)
    if body is not None:
        c.putheader("Content-Length", str(len(body)))
    c.endheaders(body)
    r = c.getresponse()
    data = r.read()
    c.close()
    return r.status, data


def test_web_rejects_foreign_hosts_tokens_and_traversal(workbench, tmp_path):
    port = workbench.server_port
    assert _req(port, "GET", "/api/session", host="evil.com")[0] == 403
    assert _req(port, "GET", "/api/session", host=f"evil.com:{port}")[0] == 403
    assert _req(port, "POST", "/api/run", b"x")[0] == 403
    for path in ("/api/files/../../etc/passwd", "/api/files/%2e%2e/x/y", "/../pyproject.toml"):
        assert _req(port, "GET", path)[0] == 404
    token = json.loads(_req(port, "GET", "/api/session")[1])["token"]
    status, _ = _req(port, "POST", "/api/run", b"hi", {"X-Workbench-Token": token,
                                                      "X-Options": quote(json.dumps({"action": "encode"})),
                                                      "X-Filename": quote("../../evil.txt")})
    assert status == 202
    assert all(p.parent.name == "input" for p in tmp_path.rglob("evil.txt"))
    assert not (tmp_path.parent / "evil.txt").exists()


# --------------------------------------------------------------------------- known bugs


@pytest.mark.xfail(strict=True, reason="BUG: a config file without `primers:` resolves the default "
                                       "primers.yaml next to the config instead of the bundled primers")
def test_partial_config_uses_bundled_primers(tmp_path):
    cfg_path = tmp_path / "my.yaml"
    cfg_path.write_text("oligo_len: 200\n")
    assert Path(load_config(cfg_path)["primers"]).exists()


@pytest.mark.xfail(strict=True, reason="BUG: non-integer oligo_len escapes the CLI as a TypeError traceback")
def test_cli_reports_bad_config_types_cleanly(tmp_path):
    (tmp_path / "primers.yaml").write_text((ROOT / "configs" / "primers.yaml").read_text())
    (tmp_path / "bad.yaml").write_text("oligo_len: abc\n")
    (tmp_path / "in.txt").write_text("hello")
    res = CliRunner().invoke(app, ["encode", str(tmp_path / "in.txt"), "--config", str(tmp_path / "bad.yaml"),
                                   "--out-dir", str(tmp_path / "o")])
    assert res.exit_code == 2 and isinstance(res.exception, SystemExit)


@pytest.mark.xfail(strict=True, reason="BUG: reedsolo is imported by the default `repair` decoder "
                                       "(ecc._ee_column) but is only a [test] extra")
def test_reedsolo_is_a_runtime_dependency():
    import tomllib

    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["dependencies"]
    assert any(d.split("=")[0].strip().lower() == "reedsolo" for d in deps)


@pytest.mark.xfail(strict=True, reason="BUG: _ee_column trusts reedsolo's global GF tables when "
                                       "gf_exp[1] == 2, so another primitive polynomial silently breaks it")
def test_ee_decoding_survives_foreign_reedsolo_tables():
    import reedsolo

    rng = random.Random(0)
    msg = np.frombuffer(rng.randbytes(40), np.uint8)
    cw = np.concatenate([msg, rs_encode_matrix(msg[:, None], 10)[:, 0]])
    word = cw.copy()
    word[[3, 17]] ^= 0x5A
    try:
        reedsolo.init_tables(prim=0x11B, generator=2, c_exp=8)
        out = _ee_column(word, [], 10)
    finally:
        reedsolo.init_tables(prim=0x11D, generator=2, c_exp=8)
    assert out is not None and bytes(out) == bytes(cw)


@pytest.mark.xfail(strict=True, reason="BUG: fountain encoder never checks decodability; with default "
                                       "overhead most small pools (k~7) fail even from perfect reads")
def test_fountain_pool_decodes_from_perfect_reads(primers):
    rng = random.Random(0)
    failures = 0
    for seed in range(8):
        data = rng.randbytes(200)
        res = encode_bytes(data, primers, make_cfg("fountain", {}, global_seed=seed))
        out = decode_reads([s for _, s in res.oligos], primers, 200, strength="erasure", grouping="index")
        failures += not (out.ok and out.data == data)
    assert failures == 0


@pytest.mark.xfail(strict=True, reason="BUG: gzip is detected by lowercase '.gz' suffix only")
def test_read_fasta_uppercase_gz_suffix(tmp_path):
    p = tmp_path / "reads.FA.GZ"
    p.write_bytes(gzip.compress(b">a\nACGT\n"))
    assert read_fasta(str(p)) == ["ACGT"]


@pytest.mark.xfail(strict=True, reason="BUG: truncated gzip raises EOFError, which the CLI does not catch")
def test_read_fasta_truncated_gzip_is_a_user_error(tmp_path):
    p = tmp_path / "reads.fa.gz"
    p.write_bytes(gzip.compress(b">a\nACGT\n" * 50)[:-6])
    with pytest.raises((ValueError, OSError)):
        read_fasta(str(p))


@pytest.mark.xfail(strict=True, reason="BUG: web accepts unbounded mean_coverage; one request can "
                                       "exhaust host memory in the read simulator")
def test_web_bounds_simulated_coverage():
    from dnastore.web import validate_options

    with pytest.raises(ValueError):
        validate_options({"mode": "oligo", "simulate": True, "channel": {"mean_coverage": 1e7}})


@pytest.mark.xfail(strict=True, reason="BUG: web allows oligo_len 100, where the header of 5 of 6 "
                                       "codecs needs more than 15 fragments")
def test_web_minimum_oligo_len_encodes_every_codec(primers):
    from dnastore.web import CODECS, validate_options

    for key, (name, params) in CODECS.items():
        opts = validate_options({"mode": "oligo", "codec": key, "oligo_len": 100})
        cfg = load_config(overrides={"oligo_len": opts["oligo_len"], "codec": {"name": name, "params": params}})
        encode_bytes(b"hello world" * 20, primers, cfg)


@pytest.mark.xfail(strict=True, reason="BUG: cluster grouping takes the first-seen header fragment, so a "
                                       "minority contaminant pool can block recovery of the majority file")
def test_cluster_grouping_prefers_majority_header(primers):
    a, b = random.Random(1).randbytes(2000), random.Random(2).randbytes(2000)
    ra = [s for _, s in encode_bytes(a, primers, make_cfg("goldman", {})).oligos]
    rb = [s for _, s in encode_bytes(b, primers, make_cfg("goldman", {})).oligos]
    out = decode_reads(ra + rb * 5, primers, 200, grouping="cluster")
    assert out.ok and out.data == b
