"""Container: file-type detection, compression decision, and the header format.

The file is treated as an opaque byte stream. Media are never decoded to pixels or
samples. The file type is detected from magic bytes (not the extension) and stored in
the header so that the decoder can restore a sensible extension.

Compression rule (DESIGN.md F12): zstd is *tried* on every file and the ratio is always
logged. The compressed form is stored only if the type is not already compressed **and**
``compressed_len / original_len < max_ratio`` (default 0.98).
"""

from __future__ import annotations

import json
import struct
import zlib
from dataclasses import dataclass, field
from typing import Any

import zstandard

from . import FORMAT_VERSION

# code: (name, extension, already_compressed)
FILE_TYPES: dict[int, tuple[str, str, bool]] = {
    0: ("bin", ".bin", False),
    1: ("txt", ".txt", False),
    2: ("png", ".png", True),
    3: ("jpeg", ".jpg", True),
    4: ("gif", ".gif", True),
    5: ("wav", ".wav", False),
    6: ("mp3", ".mp3", True),
    7: ("mp4", ".mp4", True),
    8: ("pdf", ".pdf", False),
    9: ("zip", ".zip", True),
    10: ("flac", ".flac", True),
    11: ("ogg", ".ogg", True),
    12: ("webp", ".webp", True),
}
_NAME_TO_CODE = {v[0]: k for k, v in FILE_TYPES.items()}


def detect_file_type(data: bytes) -> int:
    """Detect a file-type code from magic bytes. Falls back to ``txt`` for printable UTF-8,
    else ``bin``. Empty input is ``bin``."""
    if not data:
        return _NAME_TO_CODE["bin"]
    if data.startswith(b"\x89PNG\r\n\x1a\n"):
        return _NAME_TO_CODE["png"]
    if data.startswith(b"\xff\xd8\xff"):
        return _NAME_TO_CODE["jpeg"]
    if data[:6] in (b"GIF87a", b"GIF89a"):
        return _NAME_TO_CODE["gif"]
    if data[:4] == b"RIFF" and data[8:12] == b"WAVE":
        return _NAME_TO_CODE["wav"]
    if data[:4] == b"RIFF" and data[8:12] == b"WEBP":
        return _NAME_TO_CODE["webp"]
    if data[4:8] == b"ftyp":
        return _NAME_TO_CODE["mp4"]
    if data.startswith(b"ID3") or (len(data) > 1 and data[0] == 0xFF and (data[1] & 0xE0) == 0xE0):
        return _NAME_TO_CODE["mp3"]
    if data.startswith(b"%PDF"):
        return _NAME_TO_CODE["pdf"]
    if data.startswith(b"PK\x03\x04"):
        return _NAME_TO_CODE["zip"]
    if data.startswith(b"fLaC"):
        return _NAME_TO_CODE["flac"]
    if data.startswith(b"OggS"):
        return _NAME_TO_CODE["ogg"]
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return _NAME_TO_CODE["bin"]
    printable = sum(c.isprintable() or c in "\r\n\t" for c in text)
    return _NAME_TO_CODE["txt"] if printable >= 0.95 * len(text) else _NAME_TO_CODE["bin"]


def extension_for(code: int) -> str:
    return FILE_TYPES.get(code, FILE_TYPES[0])[1]


@dataclass
class ContainerInfo:
    """What the container stage did to the file (all of it is logged)."""

    file_type: int
    original_len: int
    stored_len: int
    compressed: bool
    trial_ratio: float  # zstd size / original size, always measured (1.0 for empty input)
    crc32: int
    reason: str


def prepare(data: bytes, zstd_level: int = 19, max_ratio: float = 0.98) -> tuple[bytes, ContainerInfo]:
    """Detect type and decide on compression. Returns ``(stored_bytes, info)``."""
    ftype = detect_file_type(data)
    crc = zlib.crc32(data) & 0xFFFFFFFF
    comp = zstandard.ZstdCompressor(level=zstd_level, write_content_size=True, write_checksum=False).compress(data)
    ratio = len(comp) / len(data) if data else 1.0
    already = FILE_TYPES[ftype][2]
    if already:
        use, reason = False, f"{FILE_TYPES[ftype][0]} is an already-compressed format"
    elif ratio < max_ratio:
        use, reason = True, f"zstd ratio {ratio:.3f} < {max_ratio}"
    else:
        use, reason = False, f"zstd ratio {ratio:.3f} >= {max_ratio}"
    stored = comp if use else data
    return stored, ContainerInfo(ftype, len(data), len(stored), use, ratio, crc, reason)


def restore(stored: bytes, compressed: bool) -> bytes:
    """Undo :func:`prepare` (decompression only; CRC is checked by the caller)."""
    if not compressed:
        return stored
    return zstandard.ZstdDecompressor().decompress(stored)


# --------------------------------------------------------------------------- header

MAGIC = b"DS"
# magic, version, file_type, flags, orig_len, stored_len, crc32, global_seed,
# codec_id, payload_bits, body_nt, seed_trits, rs_k_max, rs_parity_permille, params_len
_FIXED = struct.Struct(">2sBBBIIIIBHHBBHH")
FLAG_COMPRESSED = 0x01


@dataclass
class Header:
    """Everything the decoder needs that is not physical (primers, oligo length)."""

    file_type: int
    original_len: int
    stored_len: int
    compressed: bool
    crc32: int
    global_seed: int
    codec_id: int
    payload_bits: int
    body_nt: int
    seed_trits: int
    rs_k_max: int
    rs_parity_permille: int
    codec_params: dict[str, Any] = field(default_factory=dict)
    version: int = FORMAT_VERSION

    def pack(self) -> bytes:
        params = json.dumps(self.codec_params, separators=(",", ":"), sort_keys=True).encode()
        flags = FLAG_COMPRESSED if self.compressed else 0
        fixed = _FIXED.pack(
            MAGIC, self.version, self.file_type, flags, self.original_len, self.stored_len,
            self.crc32, self.global_seed, self.codec_id, self.payload_bits, self.body_nt,
            self.seed_trits, self.rs_k_max, self.rs_parity_permille, len(params),
        )
        return fixed + params

    @classmethod
    def unpack(cls, blob: bytes) -> Header:
        """Parse a header; trailing padding is ignored. Raises ``ValueError`` if malformed."""
        if len(blob) < _FIXED.size:
            raise ValueError("header too short")
        (magic, version, ftype, flags, olen, slen, crc, seed, codec, pb, body, strits, kmax, perm, plen) = _FIXED.unpack_from(blob)
        if magic != MAGIC:
            raise ValueError("bad magic")
        if version != FORMAT_VERSION:
            raise ValueError(f"unsupported format version {version}")
        raw = blob[_FIXED.size : _FIXED.size + plen]
        if len(raw) != plen:
            raise ValueError("header truncated")
        params = json.loads(raw.decode()) if plen else {}
        if not isinstance(params, dict):
            raise ValueError("codec parameters must be a JSON object")
        if flags & ~FLAG_COMPRESSED:
            raise ValueError("unsupported header flags")
        return cls(ftype, olen, slen, bool(flags & FLAG_COMPRESSED), crc, seed, codec, pb, body, strits, kmax, perm,
                   params, version)
