"""Regression checks for file IO, invalid configurations and the public CLI."""
import gzip
import json
import random

import pytest
from typer.testing import CliRunner

from dnastore.channel import ChannelParams, simulate
from dnastore.cli import app
from dnastore.config import load_config
from dnastore.container import Header
from dnastore.decoder import decode_reads, read_fasta
from dnastore.ecc import plan_layout
from dnastore.encoder import encode_bytes


def test_wrapped_fastq_gzip_and_fasta(tmp_path):
    fq = tmp_path / "reads.fastq.gz"
    with gzip.open(fq, "wt") as fh:
        fh.write("\n@r1\nacgt\nNN\n+\n@@@@\n++\n@r2\nTGCA\n+\n!!!!\n")
    assert read_fasta(str(fq)) == ["ACGTNN", "TGCA"]
    fa = tmp_path / "pool.fa"
    fa.write_text("\n>one\nacgt\nTG\n\n>two\nCA\n")
    assert read_fasta(str(fa)) == ["ACGTTG", "CA"]


@pytest.mark.parametrize("text", ["ACGT", ">empty\n", "@r\nACGT\n+\n!!\n", "@r\nACGT\n+\n!!!!!\n",
                                  "@r\nACGT\n", ">r\nXYZ\n", ">a\n>b\nACGT\n"])
def test_malformed_reads_rejected(tmp_path, text):
    path = tmp_path / "bad.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_fasta(str(path))


def test_ambiguous_reads_are_counted(primers):
    cfg = load_config()
    data = random.Random(0).randbytes(80)
    enc = encode_bytes(data, primers, cfg)
    result = decode_reads(["N" * 200] + [s.lower() for _, s in enc.oligos], primers, 200)
    assert result.ok and result.data == data
    assert result.report["ambiguous_reads"] == 1


@pytest.mark.parametrize("kwargs", [{"p_sub": -0.1}, {"p_ins": float("nan")}, {"mean_coverage": -1},
                                    {"p_sub": 0.8, "p_del": 0.3}, {"dropout": 2}])
def test_invalid_channel_rejected(kwargs):
    with pytest.raises(ValueError):
        ChannelParams(**kwargs)


@pytest.mark.parametrize("stored,bits,parity,kmax", [(1, 20, -1, None), (-1, 20, 150, None),
                                                    (1, 20, 150, 0), (1, 20, 150, -1)])
def test_invalid_layout_rejected(stored, bits, parity, kmax):
    with pytest.raises(ValueError):
        plan_layout(stored, bits, parity, kmax)


def test_zero_header_copies_rejected(primers):
    with pytest.raises(ValueError, match="header_copies"):
        encode_bytes(b"hello", primers, load_config(overrides={"header_copies": 0}))


def test_header_params_must_be_mapping():
    header = Header(0, 1, 1, False, 0, 0, 1, 200, 160, 3, 200, 150, [])
    with pytest.raises(ValueError, match="JSON object"):
        Header.unpack(header.pack())


def test_cli_roundtrip_outside_checkout_and_strict_baseline(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    data = random.Random(100).randbytes(400)
    path = tmp_path / "input.bin"
    path.write_bytes(data)
    runner = CliRunner()
    args = ["encode", str(path), "--codec", "naive2bit", "--out-dir", "encoded"]
    strict = runner.invoke(app, args)
    assert strict.exit_code == 2, strict.output
    assert not list((tmp_path / "encoded").glob("*.fasta"))
    manifest = json.loads(next((tmp_path / "encoded").glob("*.json")).read_text())
    assert manifest["screening"]["n_violating_oligos"] > 0
    encoded = runner.invoke(app, args + ["--allow-violations"])
    assert encoded.exit_code == 0, encoded.output
    fasta = next((tmp_path / "encoded").glob("*.fasta"))
    decoded = runner.invoke(app, ["decode", str(fasta), "--out", "recovered", "--strength", "repair"])
    assert decoded.exit_code == 0, decoded.output
    assert (tmp_path / "recovered.bin").read_bytes() == data


def test_cli_input_errors_are_concise(tmp_path):
    path = tmp_path / "bad.fastq"
    path.write_text("@r\nACGT\n+\n!!\n")
    result = CliRunner().invoke(app, ["decode", str(path)])
    assert result.exit_code == 2 and "truncated" in result.output
    assert "Traceback" not in result.output


@pytest.mark.parametrize("name,params", [("steering", {"P": -1}), ("steering", {"beam_width": 0}),
                                        ("fountain", {"overhead": -1}), ("fountain", {"delta": 1}),
                                        ("naive2bit", {"typo": True})])
def test_invalid_codec_settings_rejected(name, params):
    from dnastore.codecs import make_codec
    with pytest.raises(ValueError):
        make_codec(name, params)


def test_cli_invalid_yaml(tmp_path):
    config = tmp_path / 'bad.yaml'
    config.write_text('codec: [\n')
    payload = tmp_path / 'x.txt'
    payload.write_text('hi')
    result = CliRunner().invoke(app, ['encode', str(payload), '--config', str(config)])
    assert result.exit_code == 2 and 'Error:' in result.output


def test_valid_header_groups_do_not_consume_data_repair_budget(primers):
    cfg = load_config()
    data = random.Random(7).randbytes(100)
    enc = encode_bytes(data, primers, cfg)
    reads = [s for _, s in enc.oligos for _ in range(3)]
    result = decode_reads(reads, primers, 200, strength='repair')
    assert result.ok and result.data == data
    assert result.report['repair_attempts'] == 0
    assert result.report['data_crc_fail'] == 0


def test_mixed_header_data_cluster_is_preserved(primers, monkeypatch):
    from dnastore import decoder
    cfg = load_config()
    data = random.Random(8).randbytes(100)
    enc = encode_bytes(data, primers, cfg)
    reads = [s for _, s in enc.oligos for _ in range(3)]
    monkeypatch.setattr(decoder, 'kmer_cluster', lambda bodies, **_: [list(range(len(bodies)))])
    result = decode_reads(reads, primers, 200, strength='repair')
    assert result.ok and result.data == data
    assert result.report['groups_multi_index'] > 0
