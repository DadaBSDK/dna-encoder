"""Baseline: Goldman-style rotating ternary, no steering (DESIGN.md D5).

The frame (``floor(payload_nt * log2 3)`` bits) is converted to exactly ``payload_nt`` trits and
written with the rotating code chained from the last index base, so the payload has no
homopolymers. This is the *code* only. Goldman et al.'s 4x overlapping segments are a
redundancy scheme and are not reproduced here: all codecs share the same outer RS.
Equivalent to steering-B with P = inf. Seed 0 is the plain code; seeds > 0 XOR the frame
with a keystream first (used only for framework rerolls).
"""

from __future__ import annotations

from typing import Any

from ..dnautil import int_to_trits, keystream_int, max_bits_for_trits, rotate_decode, rotate_encode, trits_to_int
from . import Codec, EncodedPayload, register

REROLL_DOMAIN = 2


def _mask(v: int, nbits: int, index: int, seed: int, global_seed: int) -> int:
    if seed == 0:
        return v
    return v ^ keystream_int(nbits, global_seed, REROLL_DOMAIN, index, seed)


@register
class Goldman(Codec):
    codec_id = 2
    name = "goldman"

    def frame_bits(self, payload_nt: int) -> int:
        return max_bits_for_trits(payload_nt)

    def encode_payload(self, frame: int, index: int, seed: int, left: str, right: str, payload_nt: int,
                       global_seed: int) -> EncodedPayload:
        v = _mask(frame, self.frame_bits(payload_nt), index, seed, global_seed)
        return EncodedPayload(rotate_encode(int_to_trits(v, payload_nt), left[-1]))

    def decode_payload(self, seq: str, index: int, seed: int, prev: str, global_seed: int) -> int | None:
        trits = rotate_decode(seq, prev)
        if trits is None:
            return None
        v, nbits = trits_to_int(trits), max_bits_for_trits(len(seq))
        return None if v >> nbits else _mask(v, nbits, index, seed, global_seed)

    def params(self) -> dict[str, Any]:
        return {}
