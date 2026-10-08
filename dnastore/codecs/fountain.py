"""Baseline: our own LT fountain code with screening (DNA-Fountain-style; not their code).

Pipeline (replaces the outer RS for this arm; everything else is shared):
stored bytes -> whitened equal-length segments -> LT droplets -> CRC-32 -> whitened 2-bit
mapping -> screening -> keep or discard.

* **Droplet identity:** the oligo index (shared 15-nt index field) determines the droplet
  seed, ``seed32 = H(global_seed, index)``. The seed field is always 0, so screening works
  as in Erlich & Zielinski (2017): a droplet whose oligo fails a hard constraint is
  *discarded* and the next index is tried. Rejections are counted and reported.
* **Screening** applies the same hard constraints as every other arm (screening config:
  homopolymer, GC global/window, primer, 3'-anchored primer).
* **Segments** are whitened before LT encoding. The LT core's notes warn that screening
  biases which droplets survive, so whitening keeps data patterns from correlating with
  rejection.
* **Redundancy** is set by ``overhead`` (droplets = ceil(k * (1 + overhead))). For
  equal-total-nt comparisons, the benchmark sets it to match the other arms' pool size.
* **Quantisation:** LT segments are whole bytes, so up to 7 frame bits per oligo are padding.
  This is counted against density.

LT core: :mod:`dnastore.lt` (robust soliton c=0.025, delta=0.001; peeling + GF(2) Gaussian
elimination).
"""

from __future__ import annotations

from typing import Any
import math

from ..dnautil import base4_to_int, int_to_base4, keystream, keystream_int, xor_bytes
from ..ecc import CRC_BITS
from . import Codec, EncodedPayload, register

DOMAIN_SEED = 20
DOMAIN_SEGMENT = 21
DOMAIN_MAP = 22


def droplet_seed(global_seed: int, index: int) -> int:
    return int.from_bytes(keystream(4, global_seed, DOMAIN_SEED, index), "big")


@register
class Fountain(Codec):
    codec_id = 3
    name = "fountain"
    rateless = True

    def __init__(self, overhead: float = 0.10, c: float = 0.025, delta: float = 0.001, k: int = 0,
                 seg_bytes: int = 0, max_index_tries: int = 50) -> None:
        self.overhead = float(overhead)
        self.c, self.delta = float(c), float(delta)
        self.k, self.seg_bytes = int(k), int(seg_bytes)
        self.max_index_tries = int(max_index_tries)  # rejected droplets allowed per accepted one
        if not math.isfinite(self.overhead) or self.overhead < 0:
            raise ValueError("fountain overhead must be finite and nonnegative")
        if not math.isfinite(self.c) or self.c <= 0 or not 0 < self.delta < 1:
            raise ValueError("fountain requires c > 0 and 0 < delta < 1")
        if self.k < 0 or self.seg_bytes < 0 or self.max_index_tries < 1:
            raise ValueError("invalid fountain segment count, size or retry limit")

    def frame_bits(self, payload_nt: int) -> int:
        return 2 * payload_nt

    def seg_bytes_for(self, payload_nt: int) -> int:
        return (self.frame_bits(payload_nt) - CRC_BITS) // 8

    def encode_payload(self, frame: int, index: int, seed: int, left: str, right: str, payload_nt: int,
                       global_seed: int) -> EncodedPayload:
        nbits = self.frame_bits(payload_nt)
        return EncodedPayload(int_to_base4(frame ^ keystream_int(nbits, global_seed, DOMAIN_MAP, index), payload_nt))

    def decode_payload(self, seq: str, index: int, seed: int, prev: str, global_seed: int) -> int | None:
        nbits = 2 * len(seq)
        return base4_to_int(seq) ^ keystream_int(nbits, global_seed, DOMAIN_MAP, index)

    def params(self) -> dict[str, Any]:
        return {"overhead": self.overhead, "c": self.c, "delta": self.delta, "k": self.k, "seg_bytes": self.seg_bytes}

    # ------------------------------------------------------------------ segments
    @staticmethod
    def segments(stored: bytes, seg_bytes: int, global_seed: int) -> list[bytes]:
        k = -(-len(stored) // seg_bytes)
        padded = stored + bytes(k * seg_bytes - len(stored))
        return [xor_bytes(padded[i * seg_bytes : (i + 1) * seg_bytes], keystream(seg_bytes, global_seed, DOMAIN_SEGMENT, i))
                for i in range(k)]

    @staticmethod
    def unsegment(segs: list[bytes], stored_len: int, global_seed: int) -> bytes:
        raw = b"".join(xor_bytes(s, keystream(len(s), global_seed, DOMAIN_SEGMENT, i)) for i, s in enumerate(segs))
        return raw[:stored_len]
