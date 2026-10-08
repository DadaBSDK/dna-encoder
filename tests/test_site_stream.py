"""The browser port (site/stream.js) must stay interchangeable with dnastore.stream."""

from __future__ import annotations

import random
import shutil
import subprocess

import pytest

from dnastore.stream import decode_file, encode_file

from .conftest import ROOT

NODE = shutil.which("node")
CLI = ROOT / "site" / "tests" / "cli.mjs"
pytestmark = pytest.mark.skipif(NODE is None, reason="Node.js is not installed")


def _js(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([NODE, str(CLI), *map(str, args)], capture_output=True, text=True, timeout=120)


@pytest.mark.parametrize("size", [0, 1, 255, 256, 257, 65535, 65536, 65537, 150_001])
def test_js_archives_are_byte_identical_and_cross_decode(tmp_path, size):
    src = tmp_path / "input.bin"
    src.write_bytes(random.Random(size).randbytes(size))
    encode_file(src, tmp_path / "py.dna")
    assert _js("encode", src, tmp_path / "js.dna", src.name).returncode == 0
    assert (tmp_path / "js.dna").read_bytes() == (tmp_path / "py.dna").read_bytes()
    assert _js("decode", tmp_path / "py.dna", tmp_path / "js.out").returncode == 0
    assert (tmp_path / "js.out").read_bytes() == src.read_bytes()
    decode_file(tmp_path / "js.dna", tmp_path / "py.out")
    assert (tmp_path / "py.out").read_bytes() == src.read_bytes()


@pytest.mark.parametrize("name", ['quote"back\\slash.txt', "tab\tnew\nline", "ünïcødé.png", "emoji 🧬.gif",
                                  "del\x7fchar", "ctrl\x01\x1f.bin", "x" * 200])
def test_js_escapes_filenames_like_python(tmp_path, name):
    src = tmp_path / name
    src.write_bytes(b"hello DNA")
    encode_file(src, tmp_path / "py.dna")
    assert _js("encode", src, tmp_path / "js.dna", name).returncode == 0
    assert (tmp_path / "js.dna").read_bytes() == (tmp_path / "py.dna").read_bytes()


def test_js_and_python_agree_on_corrupted_archives(tmp_path):
    rng = random.Random(0)
    src = tmp_path / "x.bin"
    src.write_bytes(rng.randbytes(3000))
    encode_file(src, tmp_path / "x.dna")
    blob = (tmp_path / "x.dna").read_bytes()
    variants = [blob[:rng.randrange(len(blob))] for _ in range(25)]
    for _ in range(50):
        b = bytearray(blob)
        b[rng.randrange(len(b))] ^= 1 << rng.randrange(8)
        variants.append(bytes(b))
    variants += [blob + b"x", blob.replace(b"\n", b"\r\n"), blob.lower(), b"", b">DNASTORE-STREAM-1 [1]\n"]
    for n, v in enumerate(variants):
        (tmp_path / "m.dna").write_bytes(v)
        try:
            decode_file(tmp_path / "m.dna", tmp_path / f"py{n}")
            py_ok = True
        except ValueError:
            py_ok = False
        js = _js("decode", tmp_path / "m.dna", tmp_path / f"js{n}")
        assert js.returncode in (0, 3), js.stderr
        assert (js.returncode == 0) == py_ok, (n, js.stderr)
