import hashlib
import json

import pytest
from typer.testing import CliRunner

from dnastore.cli import app
from dnastore.stream import decode_file, encode_file


@pytest.mark.parametrize('data', [b'', bytes(range(256)), bytes(range(256)) * 1025])
def test_exact_stream_roundtrip(tmp_path, data):
    original = tmp_path / 'original.gif'
    archive, recovered = tmp_path / 'archive.dna', tmp_path / 'recovered.gif'
    original.write_bytes(data)
    enc = encode_file(original, archive)
    assert decode_file(archive, recovered) == enc
    assert recovered.read_bytes() == data
    assert enc['nucleotides'] == 4 * len(data)
    assert enc['sha256'] == hashlib.sha256(data).hexdigest()


@pytest.mark.parametrize('damage', ['substitution', 'truncated', 'tail', 'invalid', 'length', 'metadata'])
def test_corruption_never_publishes_output(tmp_path, damage):
    original, archive, recovered = (tmp_path / name for name in ['in', 'archive', 'out'])
    original.write_bytes(bytes(range(256)))
    encode_file(original, archive)
    lines = archive.read_bytes().splitlines(keepends=True)
    if damage == 'substitution':
        lines[1] = b'T' + lines[1][1:]
    elif damage == 'truncated':
        lines.pop()
    elif damage == 'tail':
        lines.append(b'A\n')
    elif damage == 'invalid':
        lines[1] = b'N' + lines[1][1:]
    elif damage == 'length':
        lines[1] = lines[1][1:]
    else:
        lines[0] = lines[0].replace(b'"in"', b'"different"')
    archive.write_bytes(b''.join(lines))
    with pytest.raises(ValueError):
        decode_file(archive, recovered)
    assert not recovered.exists()
    assert not list(tmp_path.glob('.dnastore-*'))


def test_existing_files_preserved(tmp_path):
    original, archive = tmp_path / 'in', tmp_path / 'archive'
    original.write_bytes(b'abc')
    encode_file(original, archive)
    before = archive.read_bytes()
    with pytest.raises(ValueError):
        encode_file(original, archive)
    with pytest.raises(ValueError):
        decode_file(archive, original)
    assert original.read_bytes() == b'abc'
    assert archive.read_bytes() == before


def test_cli_media_extension_and_paths(tmp_path):
    source, archive, dest = (tmp_path / s for s in ['animation.gif', 'archive.dna', 'restored.gif'])
    source.write_bytes(b'GIF89a' + bytes(range(256)))
    runner = CliRunner()
    encoded = runner.invoke(app, ['stream-encode', str(source), '--out', str(archive)])
    assert encoded.exit_code == 0, encoded.output
    decoded = runner.invoke(app, ['stream-decode', str(archive), '--out', str(dest)])
    assert decoded.exit_code == 0, decoded.output
    assert json.loads(decoded.output)['name'] == 'animation.gif'
    assert dest.read_bytes() == source.read_bytes()
