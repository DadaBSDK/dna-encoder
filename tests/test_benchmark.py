import json
import random

import pytest
import yaml

from dnastore.benchmark import depth_boundary, derived_seed, match_budget, run_benchmark, wilson
from dnastore.config import load_config
from dnastore.decoder import decode_reads
from dnastore.encoder import encode_bytes

FAST = {"w_conf": 0, "beam_width": 2, "rescore": None, "normalize": False}
ARMS = [
    {"name": "whiten", "codec": {"name": "naive2bit", "params": {"whiten": True}}},
    {"name": "goldman", "codec": {"name": "goldman"}},
    *[{"name": a, "codec": {"name": "steering", "params": {"P": 8, "alphabet": a, **FAST}}} for a in "ABC"],
    {"name": "fountain", "codec": {"name": "fountain", "params": {"overhead": .15}}},
]


def test_exact_budget_and_roundtrip_all_arms(primers):
    data = random.Random(5).randbytes(800)
    cfg = load_config()
    n, configs = match_budget(data, primers, cfg, ARMS)
    results = [encode_bytes(data, primers, configs[a["name"]]) for a in ARMS]
    assert {sum(len(s) for _, s in r.oligos) for r in results} == {n * cfg["oligo_len"]}
    for result in results:
        decoded = decode_reads([s for _, s in result.oligos], primers, 200)
        assert decoded.ok and decoded.data == data, (result.codec.name, decoded.report)
    assert configs["whiten"]["ecc"]["parity_permille"] > configs["B"]["ecc"]["parity_permille"]


def test_seed_stability_and_wilson():
    assert derived_seed(4, "arm", 3) == derived_seed(4, "arm", 3)
    assert derived_seed(4, "arm", 3) != derived_seed(4, "arm", 4)
    assert wilson(0, 20)[0] == 0
    lo, hi = wilson(20, 20)
    assert .83 < lo < .85 and hi == 1
    assert wilson(10, 20) == pytest.approx((.299298, .700702), abs=1e-6)


def test_depth_bisection_and_censoring():
    checked = []
    def probe(d):
        checked.append(d)
        return d >= 7
    assert depth_boundary(probe, 20) == 7
    assert 6 in checked and 7 in checked and len(checked) <= 6
    assert depth_boundary(lambda _: False, 20) is None
    assert depth_boundary(lambda _: True, 20) == 1


def spec():
    return {"mode": "smoke", "seed": 12, "payloads": [{"name": "random", "bytes": 80}],
            "payload_seeds": [0], "arms": ARMS[:2], "channels": [{"name": "clean", "kind": "ids"}],
            "depths": [10], "max_depth": 10, "minimum_depth": False, "channel_seeds": 1, "plots": False}


def test_run_artifacts_reproducible_and_refuses_overwrite(tmp_path):
    config = tmp_path / "bench.yaml"
    config.write_text(yaml.safe_dump(spec()))
    out = tmp_path / "run"
    rows = run_benchmark(config, out, progress=lambda _: None)
    assert len(rows) == 4 and all(r["ok"] for r in rows)
    assert len({r["total_nt"] for r in rows}) == 1
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["status"] == "complete" and prov["source_sha256"]
    assert (out / "source.zip").exists() and prov["source_archive_sha256"]
    assert (out / "summary.csv").exists() and (out / "encode.csv").exists()
    with pytest.raises(ValueError, match="not empty"):
        run_benchmark(config, out)
    rows2 = run_benchmark(config, tmp_path / "rerun", progress=lambda _: None)
    assert [{k: v for k, v in r.items() if k != "decode_s"} for r in rows] == [
        {k: v for k, v in r.items() if k != "decode_s"} for r in rows2]


def test_study_requires_twenty_seeds(tmp_path):
    config = tmp_path / "bench.yaml"
    config.write_text(yaml.safe_dump({**spec(), "mode": "study"}))
    with pytest.raises(ValueError, match=">=20"):
        run_benchmark(config, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_benchmark_failure_keeps_provenance(tmp_path):
    config = tmp_path / 'empty.yaml'
    empty = tmp_path / 'empty.bin'
    empty.write_bytes(b'')
    config.write_text(yaml.safe_dump({**spec(), 'payloads': [{'name': 'empty', 'path': 'empty.bin'}]}))
    out = tmp_path / 'failed'
    with pytest.raises(ValueError, match='nonempty'):
        run_benchmark(config, out, progress=lambda _: None)
    provenance = json.loads((out / 'provenance.json').read_text())
    assert provenance['status'] == 'failed' and 'nonempty' in provenance['error']


def test_benchmark_depth_search_exports_both_strengths(tmp_path):
    config = tmp_path / 'depth.yaml'
    config.write_text(yaml.safe_dump({**spec(), 'minimum_depth': True, 'max_depth': 10}))
    out = tmp_path / 'search'
    rows = run_benchmark(config, out, progress=lambda _: None)
    assert {r['strength'] for r in rows} == {'erasure', 'repair'}
    assert len({r['depth'] for r in rows}) > 1
    import csv
    boundaries = list(csv.DictReader((out / 'minimum_depth.csv').open()))
    assert len(boundaries) == 4 and all(r['bisection_depth'] for r in boundaries)
