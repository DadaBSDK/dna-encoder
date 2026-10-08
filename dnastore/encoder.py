"""Encode pipeline: file bytes -> container -> RS -> CRC -> codec -> oligos.

Oligo order in the output: header oligos (index < HEADER_INDEX_SPACE), then data and
parity oligos in index order. Order carries no information; the decoder ignores it.
"""

from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Any

from .codecs import Codec, make_codec
from .container import ContainerInfo, Header, prepare
from .ecc import CRC_BITS, RSLayout, add_crc_bits, plan_layout, rs_encode
from .oligo import (
    HEADER_INDEX_SPACE,
    HEADER_SEED_TRITS,
    INDEX_TRITS,
    OligoGeometry,
    assemble,
    encode_header_payload,
    encode_prefix,
    header_fragments,
    header_index,
)
from .primers import PrimerPair
from .screening import ScreeningConfig, primer_hits, screen_oligo

# Framework-level hard constraints, enforced for EVERY arm by seed reroll.
FRAMEWORK_CHECKS = ("primer_anchor",)


@dataclass(frozen=True)
class FountainLayout:
    """Layout summary for the rateless arm (duck-types the RSLayout fields used by metrics)."""

    stored_len: int
    payload_bits: int
    k: int
    n_oligos: int
    rejected: int
    blocks: tuple = ()

    @property
    def row_bytes(self) -> int:
        return self.payload_bits // 8

    @property
    def n_data_rows(self) -> int:
        return self.k

    @property
    def n_parity_rows(self) -> int:
        return self.n_oligos - self.k

    @property
    def data_oligo_equiv(self) -> float:
        return float(self.k)


def _encode_fountain(codec, stored: bytes, primers: PrimerPair, geom: OligoGeometry, gseed: int,
                     scfg: dict[str, Any]) -> tuple[list[tuple[int, str]], list[dict[str, Any]], FountainLayout]:
    """LT droplets with DNA-Fountain-style screening: failing droplets are discarded."""
    import math

    from .codecs.fountain import Fountain, droplet_seed
    from .lt import RobustSoliton, encode_droplet

    seg_bytes = codec.seg_bytes_for(geom.payload_nt)
    segs = Fountain.segments(stored, seg_bytes, gseed) if stored else []
    k = len(segs)
    codec.k, codec.seg_bytes = k, seg_bytes
    n_needed = math.ceil(k * (1 + codec.overhead)) if k else 0
    sc = ScreeningConfig.from_dict(scfg)
    dist = RobustSoliton(k, codec.c, codec.delta) if k else None
    oligos, info, idx, rejected = [], [], HEADER_INDEX_SPACE, 0
    while len(oligos) < n_needed:
        if idx >= 3**INDEX_TRITS:
            raise ValueError("file too large for the index field")
        if rejected > codec.max_index_tries * max(1, n_needed):
            raise RuntimeError(f"fountain screening rejected {rejected} droplets; constraints too strict")
        data = encode_droplet(segs, droplet_seed(gseed, idx), dist)
        frame = add_crc_bits(idx, int.from_bytes(data, "big"), 8 * seg_bytes)
        prefix = encode_prefix(idx, 0, primers.forward[-1], geom.seed_trits)
        enc = codec.encode_payload(frame, idx, 0, primers.forward + prefix, primers.tail, geom.payload_nt, gseed)
        oligo = assemble(primers.forward, prefix + enc.seq, primers.tail)
        viol = screen_oligo(oligo, primers, sc, False)["violations"]
        if viol:
            rejected += 1
        else:
            oligos.append((idx, oligo))
            info.append({"index": idx, "seed": 0, "tries": 1, "violations": [], "seed_exhausted": False})
        idx += 1
    return oligos, info, FountainLayout(len(stored), 8 * seg_bytes, k, len(oligos), rejected)


@dataclass
class EncodeResult:
    oligos: list[tuple[int, str]]
    header: Header
    info: ContainerInfo
    layout: RSLayout
    geometry: OligoGeometry
    codec: Codec
    n_header_oligos: int
    payload_info: list[dict[str, Any]] = field(default_factory=list)
    timings: dict[str, float] = field(default_factory=dict)


def framework_violations(oligo: str, primers: PrimerPair, scfg: dict[str, Any]) -> list[str]:
    """Hard constraints that apply to every arm (currently: 3'-anchored primer matches)."""
    ph = primer_hits(oligo, primers, scfg.get("primer_max_mismatch", 3), scfg.get("primer_anchor_len", 8),
                     scfg.get("primer_anchor_mismatch", 0))
    return ["primer_anchor"] if ph["anchor"] else []


def _with_seeds(make, primers: PrimerPair, scfg: dict[str, Any], n_seeds: int,
                full_screen: bool = False) -> tuple[str, dict[str, Any]]:
    """Try seeds 0..n_seeds-1 until ``make(seed) -> (oligo, info)`` has no violations.

    ``full_screen=True`` (header oligos, which are shared infrastructure in every arm) rerolls
    against all hard constraints in the screening config, not only the framework ones.

    If every seed fails, return the candidate with the fewest violations (lowest seed on ties),
    with its violations recorded and ``seed_exhausted=True``.
    """
    best = None
    for seed in range(n_seeds):
        oligo, info = make(seed)
        viol = list(info.get("violations", [])) + framework_violations(oligo, primers, scfg)
        if full_screen:
            viol += [v for v in screen_oligo(oligo, primers, ScreeningConfig.from_dict(scfg), False)["violations"] if v not in viol]
        if best is None or len(viol) < len(best[2]):
            best = (oligo, info, viol, seed)
        if not viol:
            break
    oligo, info, viol, seed = best
    return oligo, {**info, "seed": seed, "tries": seed + 1 if not viol else n_seeds, "violations": viol,
                   "seed_exhausted": bool(viol)}


def _encode_one(args: tuple) -> tuple[int, str, dict[str, Any]]:
    codec, idx, payload, primers, geom, gseed, scfg = args
    frame = add_crc_bits(idx, payload, codec.frame_bits(geom.payload_nt) - CRC_BITS)

    def make(seed: int) -> tuple[str, dict[str, Any]]:
        prefix = encode_prefix(idx, seed, primers.forward[-1], geom.seed_trits)
        enc = codec.encode_payload(frame, idx, seed, primers.forward + prefix, primers.tail, geom.payload_nt, gseed)
        if len(enc.seq) != geom.payload_nt:
            raise RuntimeError(f"codec {codec.name} produced {len(enc.seq)} nt, expected {geom.payload_nt}")
        return assemble(primers.forward, prefix + enc.seq, primers.tail), enc.info

    oligo, info = _with_seeds(make, primers, scfg, geom.n_seeds)
    return idx, oligo, info


def encode_bytes(data: bytes, primers: PrimerPair, cfg: dict[str, Any]) -> EncodeResult:
    """Encode ``data`` into an oligo pool according to ``cfg`` (see :mod:`dnastore.config`)."""
    t0 = time.perf_counter()
    geom = OligoGeometry(cfg["oligo_len"], len(primers.forward), len(primers.reverse), int(cfg.get("seed_trits", 3)))
    geom.validate()
    copies = cfg["header_copies"]
    if not isinstance(copies, int) or not 1 <= copies <= HEADER_INDEX_SPACE:
        raise ValueError("header_copies must be an integer in [1, 64]")
    codec = make_codec(cfg["codec"]["name"], cfg["codec"].get("params"))
    codec.bind(primers, cfg.get("screening", {}))
    payload_bits = codec.frame_bits(geom.payload_nt) - CRC_BITS
    if not 8 <= payload_bits < 2**16:
        raise ValueError(f"payload bits per oligo = {payload_bits}; adjust oligo_len")
    seed = int(cfg["global_seed"])
    if not 0 <= seed < 2**32:
        raise ValueError("global_seed must fit an unsigned 32-bit integer")

    stored, info = prepare(data, cfg["container"]["zstd_level"], cfg["container"]["max_ratio"])
    t1 = time.perf_counter()
    rateless = getattr(codec, "rateless", False)
    if rateless:
        f_oligos, f_info, layout = _encode_fountain(codec, stored, primers, geom, seed, cfg.get("screening", {}))
        payload_bits = layout.payload_bits
    else:
        layout = plan_layout(len(stored), payload_bits, cfg["ecc"]["parity_permille"], cfg["ecc"].get("k_max"))
    if not rateless and HEADER_INDEX_SPACE + layout.n_oligos > 3**INDEX_TRITS:
        raise ValueError("file too large for the index field")
    header = Header(
        file_type=info.file_type, original_len=info.original_len, stored_len=info.stored_len,
        compressed=info.compressed, crc32=info.crc32, global_seed=seed, codec_id=codec.codec_id,
        payload_bits=payload_bits, body_nt=geom.body_nt, seed_trits=geom.seed_trits,
        rs_k_max=0 if rateless else layout.k_max, rs_parity_permille=0 if rateless else layout.parity_permille,
        codec_params=codec.params(),
    )

    scfg = cfg.get("screening", {})
    oligos: list[tuple[int, str]] = []
    payload_info: list[dict[str, Any]] = []
    copies = int(cfg["header_copies"])
    for f, frag in enumerate(header_fragments(header.pack(), geom.header_payload_nt)):
        for c in range(copies):
            idx = header_index(f, c, copies)

            def make(seed: int, idx=idx, frag=frag) -> tuple[str, dict[str, Any]]:
                prefix = encode_prefix(idx, seed, primers.forward[-1], HEADER_SEED_TRITS)
                body = prefix + encode_header_payload(frag, idx, seed, prefix[-1], geom.header_payload_nt)
                return assemble(primers.forward, body, primers.tail), {}

            oligo, hinfo = _with_seeds(make, primers, scfg, 3**HEADER_SEED_TRITS, full_screen=True)
            oligos.append((idx, oligo))
            payload_info.append({"index": idx, "header": True, **hinfo})
    n_header = len(oligos)

    if rateless:
        oligos.extend(f_oligos)
        payload_info.extend(f_info)
        t2 = t3 = time.perf_counter()
        return EncodeResult(oligos, header, info, layout, geom, codec, n_header, payload_info,
                            {"container_s": t1 - t0, "ecc_s": 0.0, "codec_s": t3 - t1, "total_s": t3 - t0})
    frames = rs_encode(stored, layout)
    t2 = time.perf_counter()
    jobs = [(codec, HEADER_INDEX_SPACE + q, pl, primers, geom, seed, scfg) for q, pl in enumerate(frames)]
    workers = cfg.get("workers")
    if workers and workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(workers) as ex:
            results = list(ex.map(_encode_one, jobs, chunksize=max(1, len(jobs) // (4 * workers))))
    else:
        results = [_encode_one(j) for j in jobs]
    for idx, seq, info_d in results:
        oligos.append((idx, seq))
        payload_info.append({"index": idx, **info_d})
    t3 = time.perf_counter()
    return EncodeResult(
        oligos, header, info, layout, geom, codec, n_header, payload_info,
        {"container_s": t1 - t0, "ecc_s": t2 - t1, "codec_s": t3 - t2, "total_s": t3 - t0},
    )


def write_fasta(oligos: list[tuple[int, str]], path: str) -> None:
    with open(path, "w") as fh:
        for idx, seq in oligos:
            fh.write(f">oligo_{idx}\n{seq}\n")
