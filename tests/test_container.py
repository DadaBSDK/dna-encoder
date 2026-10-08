from __future__ import annotations

import pytest

from dnastore.container import FILE_TYPES, Header, detect_file_type, prepare, restore

SAMPLES = {
    "png": b"\x89PNG\r\n\x1a\n" + bytes(30),
    "jpeg": b"\xff\xd8\xff\xe0" + bytes(30),
    "gif": b"GIF89a" + bytes(30),
    "wav": b"RIFF\x00\x00\x00\x00WAVEfmt " + bytes(30),
    "webp": b"RIFF\x00\x00\x00\x00WEBPVP8 " + bytes(30),
    "mp4": b"\x00\x00\x00\x18ftypisom" + bytes(30),
    "mp3": b"ID3\x04" + bytes(30),
    "pdf": b"%PDF-1.7\n",
    "zip": b"PK\x03\x04" + bytes(30),
    "flac": b"fLaC" + bytes(30),
    "ogg": b"OggS" + bytes(30),
    "txt": "hello, wörld\nline two\n".encode(),
    "bin": bytes(range(256)),
}


@pytest.mark.parametrize("name,data", SAMPLES.items())
def test_detect(name, data):
    assert FILE_TYPES[detect_file_type(data)][0] == name


def test_detect_ignores_extension_semantics():
    # MP3 frame sync without ID3 tag
    assert FILE_TYPES[detect_file_type(b"\xff\xfb\x90\x00" + bytes(10))][0] == "mp3"
    assert FILE_TYPES[detect_file_type(b"")][0] == "bin"


def test_compression_decision():
    text = b"the quick brown fox " * 500
    stored, info = prepare(text)
    assert info.compressed and len(stored) < len(text) and restore(stored, True) == text
    png = SAMPLES["png"] * 100  # compressible bytes, but PNG is skipped by rule
    stored, info = prepare(png)
    assert not info.compressed and stored == png and info.trial_ratio < 1.0
    import os

    rnd = os.urandom(4000)
    stored, info = prepare(rnd)
    assert not info.compressed and stored == rnd


def test_header_roundtrip():
    h = Header(3, 123456, 120000, True, 0xDEADBEEF, 42, 4, 251, 160, 4, 221, 150, {"P": 8, "alphabet": "B"})
    blob = h.pack()
    assert Header.unpack(blob + bytes(7)) == h  # trailing padding ignored


def test_header_rejects_garbage():
    with pytest.raises(ValueError):
        Header.unpack(b"XX" + bytes(40))
    with pytest.raises(ValueError):
        Header.unpack(b"DS")
