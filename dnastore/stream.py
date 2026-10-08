"""Uncapped, bounded-memory DNA text archives for exact digital file storage.

This separate format has no synthesis constraints or error correction. DNA encodes
opaque bytes with 00=A, 01=C, 10=G, 11=T; framing carries name and SHA-256.
"""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import tempfile

MAGIC = b'>DNASTORE-STREAM-1 '
BASES = b'ACGT'
ENCODE = tuple(bytes(BASES[(v >> shift) & 3] for shift in (6, 4, 2, 0)) for v in range(256))
DECODE = {seq: v for v, seq in enumerate(ENCODE)}
CHUNK_BYTES = 64 * 1024
LINE_BYTES = 256
MAX_HEADER = 16384


@contextmanager
def _new_output(path: Path):
    """Publish complete output atomically without replacing any existing file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ValueError(f'output already exists: {path}')
    fd, name = tempfile.mkstemp(prefix='.dnastore-', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as fh:
            yield fh
        os.link(name, path)  # exclusive creation, including a race with another writer
    finally:
        os.unlink(name)


def encode_file(source: str | Path, output: str | Path) -> dict:
    """Stream any file into a DNA archive. No whole-file buffering or length cap."""
    source, output = Path(source), Path(output)
    metadata = {'name': source.name, 'mapping': '00=A,01=C,10=G,11=T'}
    header = MAGIC + json.dumps(metadata, ensure_ascii=True).encode() + b'\n'
    if len(header) > MAX_HEADER:
        raise ValueError('filename metadata too long')
    digest, size = hashlib.sha256(), 0
    with source.open('rb') as src, _new_output(output) as dst:
        dst.write(header)
        while chunk := src.read(CHUNK_BYTES):
            digest.update(chunk)
            size += len(chunk)
            for start in range(0, len(chunk), LINE_BYTES):
                dst.write(b''.join(ENCODE[v] for v in chunk[start:start + LINE_BYTES]) + b'\n')
        result = {**metadata, 'bytes': size, 'nucleotides': size * 4, 'sha256': digest.hexdigest()}
        dst.write(b'!' + json.dumps(result, sort_keys=True).encode() + b'\n')
    return result


def decode_file(source: str | Path, output: str | Path) -> dict:
    """Verify framing, byte count and SHA-256 before publishing the recovered file.

    The caller chooses the output path; archive filenames are never used as paths.
    """
    source, output = Path(source), Path(output)
    with source.open('rb') as src:
        header = src.readline(MAX_HEADER + 1)
        if len(header) > MAX_HEADER or not header.endswith(b'\n') or not header.startswith(MAGIC):
            raise ValueError('not a supported DNASTORE streaming archive')
        metadata = json.loads(header[len(MAGIC):])
        if (not isinstance(metadata, dict) or not isinstance(metadata.get('name'), str)
                or metadata.get('mapping') != '00=A,01=C,10=G,11=T'):
            raise ValueError('invalid streaming metadata')
        digest, size = hashlib.sha256(), 0
        with _new_output(output) as dst:
            while True:
                line = src.readline(MAX_HEADER + 1)
                if not line:
                    raise ValueError('truncated archive: missing integrity footer')
                if len(line) > MAX_HEADER or not line.endswith(b'\n'):
                    raise ValueError('invalid or oversized archive record')
                if line.startswith(b'!'):
                    footer = json.loads(line[1:])
                    expected = {**metadata, 'bytes': size, 'nucleotides': size * 4, 'sha256': digest.hexdigest()}
                    if footer != expected:
                        raise ValueError('archive integrity check failed (metadata, length or SHA-256)')
                    if src.read(1):
                        raise ValueError('unexpected data after integrity footer')
                    return expected
                dna = line[:-1]
                if not dna or len(dna) > LINE_BYTES * 4 or len(dna) % 4:
                    raise ValueError('invalid DNA record length')
                try:
                    chunk = bytes(DECODE[dna[i:i + 4]] for i in range(0, len(dna), 4))
                except KeyError as exc:
                    raise ValueError('DNA records must contain only uppercase A, C, G, T') from exc
                dst.write(chunk)
                digest.update(chunk)
                size += len(chunk)
