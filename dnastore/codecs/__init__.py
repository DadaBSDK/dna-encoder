"""Payload codecs.

A codec maps a *frame* (an int of ``frame_bits`` bits: the oligo's share of the bit-packed
RS stream followed by the CRC-32, produced by the ECC layer) to the
oligo's payload region, i.e. the bases between the shared index+seed prefix and the
reverse primer, and back. The per-oligo ``seed`` (0 .. N_SEEDS-1) selects a re-encoding;
the encoder rerolls it when hard constraints fail. Codecs see the full left flank (primer + index) and right flank
(reverse-primer tail) so that constraint-aware codecs can score junctions.

Every codec in the RS-framed pipeline shares the primers, the index field, the inner
CRC, and the outer RS (DESIGN.md F11). Only the frame <-> bases mapping differs.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, ClassVar


@dataclass
class EncodedPayload:
    seq: str
    info: dict[str, Any] = field(default_factory=dict)


class Codec(ABC):
    """Base class for frame <-> payload-region codecs."""

    codec_id: ClassVar[int]
    name: ClassVar[str]

    @abstractmethod
    def frame_bits(self, payload_nt: int) -> int:
        """Number of frame bits (payload + CRC) carried by ``payload_nt`` bases."""

    @abstractmethod
    def encode_payload(self, frame: int, index: int, seed: int, left: str, right: str, payload_nt: int,
                       global_seed: int) -> EncodedPayload:
        """Encode a frame into exactly ``payload_nt`` bases.

        Always returns a (best-effort) sequence. If the codec's own hard constraints cannot be
        met with this seed, their names go in ``info["violations"]``. The encoder then tries
        the next seed and, if every seed fails, emits the least-violating candidate with its
        violations recorded (never silently)."""

    @abstractmethod
    def decode_payload(self, seq: str, index: int, seed: int, prev: str, global_seed: int) -> int | None:
        """Decode the payload region (``prev`` = base preceding it). ``None`` if impossible."""

    def bind(self, primers: Any, screening: dict[str, Any]) -> None:
        """Give constraint-aware codecs the primers and screening config before encoding.

        Baselines ignore this. The screening config is the single source of truth for
        constraint values (GC limits, run length, primer-match rules)."""

    @abstractmethod
    def params(self) -> dict[str, Any]:
        """Parameters stored in the header (JSON-serialisable)."""

    @classmethod
    def from_params(cls, params: dict[str, Any]) -> Codec:
        return cls(**params)


REGISTRY: dict[int, type[Codec]] = {}
BY_NAME: dict[str, type[Codec]] = {}


def register(cls: type[Codec]) -> type[Codec]:
    REGISTRY[cls.codec_id] = cls
    BY_NAME[cls.name] = cls
    return cls


def make_codec(name: str, params: dict[str, Any] | None = None) -> Codec:
    if name not in BY_NAME:
        raise ValueError(f"unknown codec {name!r}; choose from {', '.join(sorted(BY_NAME))}")
    try:
        return BY_NAME[name].from_params(params or {})
    except TypeError as exc:
        raise ValueError(f"invalid {name} codec parameters: {exc}") from exc


def codec_from_header(codec_id: int, params: dict[str, Any]) -> Codec:
    if codec_id not in REGISTRY:
        raise ValueError(f"unknown codec id {codec_id}")
    return REGISTRY[codec_id].from_params(params)


from . import fountain, goldman, naive2bit, steering  # noqa: E402,F401  (registration side effects)
