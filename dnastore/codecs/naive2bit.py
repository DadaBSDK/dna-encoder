"""Baseline: direct 2-bit mapping (00->A, 01->C, 10->G, 11->T) with no constraints.

With ``whiten=True`` the frame is XORed with a keystream keyed by (global seed, index,
seed) before mapping. That gives the unconstrained *whitened* quaternary payload used by the
primary ``whiten+RS`` comparison arm (DESIGN.md Part 0.2). The codec does no
screening of its own: constraint violations are measured, never fixed. The only exception is the
framework-level hard constraint (3'-anchored primer matches), which is enforced for every
arm by rerolling the per-oligo seed. With ``whiten=False``, seed 0 is the raw mapping and
seeds > 0 apply a keystream (needed only for rerolls).

Frames are bit strings of exactly ``2 * payload_nt`` bits, so there is no padding.
"""

from __future__ import annotations

from typing import Any

from ..dnautil import base4_to_int, int_to_base4, keystream_int
from . import Codec, EncodedPayload, register

WHITEN_DOMAIN = 1


@register
class Naive2Bit(Codec):
    codec_id = 1
    name = "naive2bit"

    def __init__(self, whiten: bool = False) -> None:
        self.whiten = bool(whiten)

    def frame_bits(self, payload_nt: int) -> int:
        return 2 * payload_nt

    def _mask(self, v: int, nbits: int, index: int, seed: int, global_seed: int) -> int:
        if not self.whiten and seed == 0:
            return v
        return v ^ keystream_int(nbits, global_seed, WHITEN_DOMAIN, index, seed)

    def encode_payload(self, frame: int, index: int, seed: int, left: str, right: str, payload_nt: int,
                       global_seed: int) -> EncodedPayload:
        nbits = self.frame_bits(payload_nt)
        return EncodedPayload(int_to_base4(self._mask(frame, nbits, index, seed, global_seed), payload_nt))

    def decode_payload(self, seq: str, index: int, seed: int, prev: str, global_seed: int) -> int | None:
        return self._mask(base4_to_int(seq), 2 * len(seq), index, seed, global_seed)

    def params(self) -> dict[str, Any]:
        return {"whiten": self.whiten}
