"""Decode pipeline: reads -> orientation/trim -> grouping by index -> consensus -> codec -> CRC -> RS -> file.

1. **Orientation and trimming** (:func:`dnastore.consensus.locate_and_orient`): edlib
   infix search for F at the read start and revcomp(R) at the end, on both strands.
2. **Grouping.** Default: similarity clustering of bodies (k-mer candidates + bounded
   edlib), independent of the index field. At nanopore error rates (~10%/nt) the 15-nt
   index decodes exactly in only ~20% of reads (Phase 3), and in rotating ternary one
   substitution corrupts two trits (DESIGN.md 0b.6). The index is then read from each
   group's candidate bodies. Alternative ``grouping="index"``: group by each read's exactly
   decoded index (fast, low-error channels only). Either way a wrong index can only cause
   a CRC failure, because the CRC is bound to the index.
3. **Per-group candidates**, tried until one passes CRC: exact bodies seen at least twice
   (most frequent first), then an alignment consensus of all bodies in the group, then a
   lone exact-length body.
4. **Fallback** (``grouping="index"`` only): reads whose index cannot be decoded are
   clustered by edit distance and each cluster is decoded like a group.
5. CRC-valid frames go to the bit-packed outer RS; then decompression and CRC32.

``DecodeResult.consensus`` maps index to the body that was decoded (or the consensus
attempted), which the benchmark compares with the true body to give the post-consensus
per-oligo error rate.
"""

from __future__ import annotations

import time
import gzip
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np

from .codecs import codec_from_header
from .consensus import cluster, consensus, kmer_cluster, locate_and_orient
from .container import Header, extension_for, restore
from .ecc import CRC_BITS, check_crc_bits, plan_layout, rs_decode
from .repair import repair
from .oligo import HEADER_INDEX_SPACE, HEADER_SEED_TRITS, INDEX_TRITS, decode_header_payload, decode_index, decode_prefix
from .primers import PrimerPair


@dataclass
class DecodeResult:
    ok: bool
    data: bytes | None
    header: Header | None
    extension: str
    report: dict[str, Any] = field(default_factory=dict)
    consensus: dict[int, str] = field(default_factory=dict)


def read_fasta(path: str) -> list[str]:
    """Read FASTA or wrapped FASTQ, optionally gzip-compressed; reject malformed records."""
    seqs: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as fh:
        first = next((line.strip() for line in fh if line.strip()), "")
        if not first:
            return []
        if first.startswith("@"):
            header = first
            while header:
                if not header.startswith("@"):
                    raise ValueError("FASTQ record must start with @")
                parts = []
                for line in fh:
                    if line.startswith("+"):
                        break
                    parts.append(line.strip())
                else:
                    raise ValueError("FASTQ record is missing its + separator")
                seq = "".join(parts).upper()
                if not seq:
                    raise ValueError("FASTQ record has an empty sequence")
                quality_len = 0
                while quality_len < len(seq):
                    line = fh.readline()
                    if not line:
                        raise ValueError("FASTQ quality is truncated")
                    quality_len += len(line.rstrip("\r\n"))
                if quality_len != len(seq):
                    raise ValueError("FASTQ quality length does not match sequence")
                seqs.append(seq)
                header = next((line.strip() for line in fh if line.strip()), "")
        elif first.startswith(">"):
            cur: list[str] = []
            for line in fh:
                if line.startswith(">"):
                    if not cur:
                        raise ValueError("FASTA record has an empty sequence")
                    seqs.append("".join(cur).upper())
                    cur = []
                elif line.strip():
                    cur.append(line.strip())
            if not cur:
                raise ValueError("FASTA record has an empty sequence")
            seqs.append("".join(cur).upper())
        else:
            raise ValueError("expected a FASTA (>) or FASTQ (@) header")
    if any(set(s) - set("ACGTN") for s in seqs):
        raise ValueError("reads must contain only A, C, G, T or N")
    return seqs


def _candidates(bodies: Counter, expected_len: int, seed: int) -> list[tuple[str, str]]:
    """Ordered (body, how) candidates for one index group."""
    out: list[tuple[str, str]] = []
    for body, n in bodies.most_common(3):
        if n >= 2 and len(body) == expected_len:
            out.append((body, "exact"))
    reads = list(bodies.elements())
    if len(reads) >= 2:
        cons, _ = consensus(reads, expected_len, rng=np.random.default_rng(seed))
        if all(cons != b for b, _ in out):
            out.append((cons, "consensus"))
    if len(reads) == 1 and len(reads[0]) == expected_len:
        out.append((reads[0], "single"))
    return out


REPAIR_READS = 50  # reads used for run-length support per group (sampled when larger)


def _group_candidates(args: tuple) -> tuple[list[tuple[str, str]], list[str]]:
    """Candidates of the whole group, then of its sub-groups by exactly decoded read index
    (so a cluster that merged several oligos can still yield each of them). Also returns
    the reads used to rank repair edits (empty for single-read groups)."""
    bodies, expected_len, seed, subgroups = args
    out = _candidates(bodies, expected_len, seed)
    seen = {b for b, _ in out}
    if len(subgroups) > 1:
        for sub in subgroups:
            for b, how in _candidates(sub, expected_len, seed):
                if b not in seen:
                    seen.add(b)
                    out.append((b, how + "/sub"))
    reads = list(bodies.elements())
    if len(reads) > REPAIR_READS:
        rng = np.random.default_rng(seed)
        reads = [reads[i] for i in rng.choice(len(reads), REPAIR_READS, replace=False)]
    return out, (reads if len(reads) >= 2 else [])


@dataclass
class DecoderOptions:
    """Decoder strength (the decoder is part of the ECC side of the comparison).

    ``strength="erasure"``: CRC-failing oligos are discarded, and RS does erasure decoding
    only (plus a search for one undetected bad oligo per block).
    ``strength="repair"``: CRC-guided homopolymer run-length repair (:mod:`dnastore.repair`)
    of failed consensus bodies, then RS errors-and-erasures decoding that also uses the
    unverified (CRC-failing) decoded payloads.
    """

    repair: bool = False
    rs_errors: bool = False
    repair_max_edits: int = 2
    repair_budget: int = 2000

    @classmethod
    def from_strength(cls, strength: str, **kw: Any) -> "DecoderOptions":
        if strength == "erasure":
            return cls(False, False, **kw)
        if strength == "repair":
            return cls(True, True, **kw)
        raise ValueError(f"unknown decoder strength {strength!r}")


# (index, payload or None, body, how, unverified payload or None, repair candidates tried, repair edits)
GroupOut = tuple[int | None, int | None, str | None, str, int | None, int, int]


def _resolve_group(args: tuple) -> list[GroupOut]:
    """One entry per distinct index that passed CRC, or a single failed entry.

    The index is read from each candidate body. ``hint`` (exact-index grouping) restricts it
    to that index; ``valid`` is the range of data indices. A wrong index cannot pass, because
    the CRC is bound to the index. On failure, with ``opts.repair`` the group's consensus is
    repaired, and with ``opts.rs_errors`` the best candidate's decoded payload is returned
    as *unverified*.
    """
    hint, (cands, reads), expected_len, codec, f_last, seed_trits, frame_bits, gseed, valid, opts = args
    pn = INDEX_TRITS + seed_trits
    nbits = frame_bits - CRC_BITS

    def parse(body: str) -> tuple[int, int] | None:
        if len(body) != expected_len:
            return None
        pre = decode_prefix(body, f_last, seed_trits)
        if pre is None or (hint is not None and pre[0] != hint) or not (valid[0] <= pre[0] < valid[1]):
            return None
        return pre

    def frame_of(body: str, pre: tuple[int, int]) -> int | None:
        fr = codec.decode_payload(body[pn:], pre[0], pre[1], body[pn - 1], gseed)
        return None if fr is None or fr < 0 or fr >> frame_bits else fr

    first = None
    found: dict[int, GroupOut] = {}
    for body, how in cands:
        pre = parse(body)
        if pre is None:
            first = first or (hint, body, None)
            continue
        first = first if first and first[2] is not None else (pre[0], body, pre)
        if pre[0] in found:
            continue
        fr = frame_of(body, pre)
        payload = None if fr is None else check_crc_bits(pre[0], fr, frame_bits)
        if payload is not None:
            found[pre[0]] = (pre[0], payload, body, how, None, 0, 0)
            if hint is not None:
                break
    if found:
        return list(found.values())
    tried = edits = 0
    if opts.repair and reads:
        base = next((b for b, h in cands if h == "consensus"), cands[0][0] if cands else None)
        hit: dict[str, Any] = {}

        def accept(body: str) -> bool:
            pre = parse(body)
            if pre is None:
                return False
            fr = frame_of(body, pre)
            payload = None if fr is None else check_crc_bits(pre[0], fr, frame_bits)
            if payload is None:
                return False
            hit.update(idx=pre[0], payload=payload)
            return True

        if base is not None:
            rr = repair(base, reads, expected_len, accept, opts.repair_max_edits, opts.repair_budget)
            tried, edits = rr.tried, rr.edits
            if rr.body is not None:
                return [(hit["idx"], hit["payload"], rr.body, "repair", None, tried, edits)]
    if first is None:
        return [(hint, None, None, "failed", None, tried, edits)]
    unverified = None
    if opts.rs_errors and first[2] is not None:
        fr = frame_of(first[1], first[2])
        unverified = None if fr is None else fr >> CRC_BITS
        if unverified is not None and unverified >> nbits:
            unverified = None
    return [(first[0], None, first[1], "failed", unverified, tried, edits)]


def _pmap(fn, jobs: list, workers: int | None) -> list:
    if workers and workers > 1 and len(jobs) > 64:
        with ProcessPoolExecutor(workers) as ex:
            return list(ex.map(fn, jobs, chunksize=max(1, len(jobs) // (8 * workers))))
    return [fn(j) for j in jobs]


def _resolve_all(groups: list[tuple[int | None, tuple]], body_nt: int, codec, f_last: str, seed_trits: int,
                 frame_bits: int, gseed: int, valid: tuple[int, int], opts: DecoderOptions,
                 workers: int | None) -> tuple[list[tuple[int | None, GroupOut]], int]:
    """Resolve every group. Also returns how many groups yielded more than one CRC-valid
    oligo (a cluster that merged oligos and was split by the per-read-index sub-groups)."""
    jobs = [(h, c, body_nt, codec, f_last, seed_trits, frame_bits, gseed, valid, opts) for h, c in groups]
    out = _pmap(_resolve_group, jobs, workers)
    multi = sum(sum(r[1] is not None for r in rs) > 1 for rs in out)
    return [(h, r) for (h, _), rs in zip(groups, out) for r in rs], multi


def decode_reads(reads: Iterable[str], primers: PrimerPair, oligo_len: int, workers: int | None = None,
                 max_primer_ed: int = 4, cluster_fallback: bool = True,
                 primer_margin: int = 80, grouping: str = "cluster",
                 strength: str | DecoderOptions = "repair") -> DecodeResult:
    """Recover the original file from a collection of reads (any order, either strand).

    ``grouping="cluster"`` (default) groups reads by similarity (:func:`kmer_cluster`) and
    reads the index from each group's candidates; it works at nanopore error rates.
    ``grouping="index"`` groups by the exactly decoded index field of each read, with
    edit-distance clustering of undecodable reads as fallback; it is only reliable when
    most reads carry an error-free index (low IDS rates).
    ``strength``: ``"erasure"`` or ``"repair"`` (see :class:`DecoderOptions`).
    """
    opts = strength if isinstance(strength, DecoderOptions) else DecoderOptions.from_strength(strength)
    t0 = time.perf_counter()
    rep: dict[str, Any] = Counter()
    f, t = primers.forward, primers.tail
    body_nt = oligo_len - len(f) - len(t)
    bodies: list[str] = []
    prov: list[int | None] = []
    for read in reads:
        rep["reads"] += 1
        read = read.upper()
        if not read or set(read) - set("ACGT"):
            rep["ambiguous_reads"] += 1
            continue
        o = locate_and_orient(read, f, t, max_primer_ed, primer_margin)
        if o is None:
            rep["unoriented"] += 1
            continue
        body, is_rev, _ = o
        rep["reverse_strand"] += int(is_rev)
        rep["exact_length"] += int(len(body) == body_nt)
        idx = decode_index(body[:INDEX_TRITS], f[-1])
        rep["bad_index"] += int(idx is None)
        bodies.append(body)
        prov.append(idx)
    rep = dict(rep)

    # ---- grouping: list of (index hint or None, Counter of bodies)
    raw_groups: list[tuple[int | None, Counter]] = []
    subs: list[list[Counter]] = []
    orphans: list[str] = []
    if grouping == "cluster":
        order = sorted(range(len(bodies)), key=lambda i: (len(bodies[i]) != body_nt, i))
        for members in kmer_cluster(bodies, order=order):
            raw_groups.append((None, Counter(bodies[m] for m in members)))
            by: dict[int, Counter] = defaultdict(Counter)
            for m in members:
                if prov[m] is not None:
                    by[prov[m]][bodies[m]] += 1
            subs.append(list(by.values()) if len(by) > 1 else [])
    elif grouping == "index":
        by_idx: dict[int, Counter] = defaultdict(Counter)
        for b, idx in zip(bodies, prov):
            if idx is None:
                orphans.append(b)
            else:
                by_idx[idx][b] += 1
        raw_groups = list(by_idx.items())
        if cluster_fallback and orphans:
            for members in cluster(orphans, rng=np.random.default_rng(0)):
                if len(members) >= 2:
                    raw_groups.append((None, Counter(orphans[m] for m in members)))
        subs = [[] for _ in raw_groups]
    else:
        raise ValueError(f"unknown grouping {grouping!r}")
    rep["groups"] = len(raw_groups)
    cands = _pmap(_group_candidates, [(g, body_nt, i, sb) for i, ((_, g), sb) in enumerate(zip(raw_groups, subs))],
                  workers)
    groups = [(h, c) for (h, _), c in zip(raw_groups, cands)]

    # ---- header (fixed format, decoded before anything codec-specific)
    frags: dict[int, bytes] = {}
    n_frags = None
    pnh = INDEX_TRITS + HEADER_SEED_TRITS

    def header_try(body: str, hint: int | None) -> tuple[int, int, bytes] | None:
        if len(body) != body_nt:
            return None
        pre = decode_prefix(body, f[-1], HEADER_SEED_TRITS)
        if pre is None or pre[0] >= HEADER_INDEX_SPACE or (hint is not None and pre[0] != hint):
            return None
        return decode_header_payload(body[pnh:], pre[0], pre[1], body[pnh - 1])

    def add_frag(got: tuple[int, int, bytes]) -> None:
        nonlocal n_frags
        fid, n, chunk = got
        n_frags = n if n_frags is None else n_frags
        frags.setdefault(fid, chunk)

    header_failed: list[tuple[int | None, tuple]] = []
    verified_header_bodies: set[str] = set()
    for hint, cl in groups:
        if hint is not None and hint >= HEADER_INDEX_SPACE:
            continue
        ok = False
        for body, _ in cl[0]:
            got = header_try(body, hint)
            if got is not None:
                verified_header_bodies.add(body)
                add_frag(got)
                ok = True
                if hint is not None:
                    break
        if not ok and cl[0]:
            base = next((b for b, h in cl[0] if h == "consensus"), cl[0][0][0])
            hi = decode_index(base[:INDEX_TRITS], f[-1])
            if (hint if hint is not None else hi) is not None and (hint if hint is not None else hi) < HEADER_INDEX_SPACE:
                header_failed.append((hint, cl))
    rep["header_repaired"] = 0
    if opts.repair and (n_frags is None or any(i not in frags for i in range(n_frags))):
        for hint, (cl, rd) in header_failed:
            if not rd:
                continue
            base = next((b for b, h in cl if h == "consensus"), cl[0][0])
            hit: list = []
            rr = repair(base, rd, body_nt, lambda b: bool(hit.append(header_try(b, hint)) or hit[-1]),
                        opts.repair_max_edits, opts.repair_budget)
            if rr.body is not None:
                add_frag(hit[-1])
                rep["header_repaired"] += 1
    rep["header_fragments"] = len(frags)
    if n_frags is None or any(i not in frags for i in range(n_frags)):
        return DecodeResult(False, None, None, ".bin", {**rep, "error": "header incomplete"})
    try:
        header = Header.unpack(b"".join(frags[i] for i in range(n_frags)))
        if header.body_nt != body_nt:
            raise ValueError("header body length does not match oligo_len and primers")
        if not 1 <= header.seed_trits <= 8 or not 8 <= header.payload_bits < 65536:
            raise ValueError("invalid seed or payload geometry")
        codec = codec_from_header(header.codec_id, header.codec_params)
        if not getattr(codec, "rateless", False):
            layout = plan_layout(header.stored_len, header.payload_bits, header.rs_parity_permille, header.rs_k_max)
        if header.payload_bits + CRC_BITS > codec.frame_bits(body_nt - INDEX_TRITS - header.seed_trits):
            raise ValueError("payload exceeds codec capacity")
    except (ValueError, TypeError, KeyError, OverflowError) as exc:
        return DecodeResult(False, None, None, ".bin", {**rep, "error": f"header parse: {exc}"})
    frame_bits = header.payload_bits + CRC_BITS
    # Cluster hints are None, including for headers. Do not spend the data-repair
    # budget trying to turn CRC-verified header-only groups into data oligos.
    # Keep mixed/uncertain groups so subgroup recovery still has all its candidates.
    data_groups = [(h, c) for h, c in groups
                   if (h is None or h >= HEADER_INDEX_SPACE)
                   and not (c[0] and all(b in verified_header_bodies for b, _ in c[0]))]
    if getattr(codec, "rateless", False):
        results, rep["groups_multi_index"] = _resolve_all(
            data_groups, body_nt, codec, f[-1], header.seed_trits, frame_bits, header.global_seed,
            (HEADER_INDEX_SPACE, 3**INDEX_TRITS), opts, workers)
        _repair_stats(rep, [r for _, r in results])
        return _decode_fountain(codec, header, [r[:4] for _, r in results], rep, t0)
    layout = plan_layout(header.stored_len, header.payload_bits, header.rs_parity_permille, header.rs_k_max)
    valid = (HEADER_INDEX_SPACE, HEADER_INDEX_SPACE + layout.n_oligos)

    # ---- data oligos
    results, rep["groups_multi_index"] = _resolve_all(data_groups, body_nt, codec, f[-1], header.seed_trits,
                                                      frame_bits, header.global_seed, valid, opts, workers)
    received: dict[int, int] = {}
    unverified: dict[int, tuple[int, int]] = {}  # position -> (group size rank key, payload)
    cons_bodies: dict[int, str] = {}
    how_counts: Counter = Counter()
    recovered_by_cluster = 0
    for gi, (hint, (idx, payload, body, how, unv, _, _)) in enumerate(results):
        how_counts[how] += 1
        if idx is None or not (valid[0] <= idx < valid[1]):
            continue
        q = idx - HEADER_INDEX_SPACE
        if payload is not None:
            if q not in received:
                received[q] = payload
                cons_bodies[idx] = body
                recovered_by_cluster += int(grouping == "index" and hint is None)
        else:
            if body is not None and q not in received:
                cons_bodies.setdefault(idx, body)
            if unv is not None and q not in unverified:
                unverified[q] = unv
    for q in received:
        unverified.pop(q, None)
    rep.update(data_oligos_expected=layout.n_oligos, data_oligos_valid=len(received),
               data_crc_fail=how_counts["failed"], resolved_by={k: v for k, v in how_counts.items()},
               recovered_by_cluster=recovered_by_cluster, indices_seen=len(cons_bodies),
               unverified_used=len(unverified) if opts.rs_errors else 0)
    _repair_stats(rep, [r for _, r in results])
    ext = extension_for(header.file_type)
    attempts = [(unverified, True)] if opts.rs_errors and unverified else []
    attempts.append(({}, False))
    last_err = ""
    for unv_map, use_errors in attempts:
        rs = rs_decode(received, layout, unv_map, errors=use_errors)
        rep.update(rs_blocks=len(layout.blocks), rs_blocks_failed=rs.blocks_failed,
                   rs_max_erasures=max(rs.erasures_per_block, default=0), rs_corrected_rows=rs.corrected_rows,
                   rs_ee_columns=rs.ee_columns, rs_errors_mode=use_errors)
        if rs.data is None:
            last_err = "RS failure"
            rep["rs_reasons"] = rs.reasons
            continue
        try:
            data = restore(rs.data, header.compressed)
        except Exception:  # zstd raises its own error type; a miscorrection can land here
            last_err = "decompress"
            continue
        crc_ok = (zlib.crc32(data) & 0xFFFFFFFF) == header.crc32 and len(data) == header.original_len
        if crc_ok:
            rep.update(crc32_ok=True, decode_s=time.perf_counter() - t0)
            return DecodeResult(True, data, header, ext, rep, cons_bodies)
        last_err = "file CRC32 mismatch"
    rep.update(crc32_ok=False, decode_s=time.perf_counter() - t0)
    return DecodeResult(False, None, header, ext, {**rep, "error": last_err}, cons_bodies)


def _repair_stats(rep: dict[str, Any], results: list[GroupOut]) -> None:
    rep["repaired"] = sum(r[3] == "repair" for r in results)
    rep["repair_attempts"] = sum(r[5] > 0 for r in results)
    rep["repair_candidates_tried"] = sum(r[5] for r in results)
    rep["repair_hit_ranks"] = [r[5] for r in results if r[3] == "repair"]
    rep["repair_hit_edits"] = [r[6] for r in results if r[3] == "repair"]


def _decode_fountain(codec, header: Header, results: list[tuple], rep: dict[str, Any], t0: float) -> DecodeResult:
    """Rateless arm: every CRC-valid droplet (any data index) goes to the LT decoder."""
    from .codecs.fountain import Fountain, droplet_seed
    from .lt import LTDecoder, RobustSoliton

    ext = extension_for(header.file_type)
    cons_bodies: dict[int, str] = {}
    for idx, payload, body, _ in results:
        if idx is not None and body is not None and (payload is not None or idx not in cons_bodies):
            cons_bodies[idx] = body
    k, sb = codec.k, codec.seg_bytes
    if k == 0:
        stored = b""
    else:
        dec = LTDecoder(k, sb, RobustSoliton(k, codec.c, codec.delta))
        n_valid = 0
        seen: set[int] = set()
        for idx, payload, _, _ in results:
            if payload is not None and idx not in seen:
                seen.add(idx)
                dec.add(droplet_seed(header.global_seed, idx), payload.to_bytes(sb, "big"))
                n_valid += 1
        res = dec.decode()
        rep.update(droplets_valid=n_valid, droplets_crc_fail=sum(p is None for _, p, _, _ in results),
                   droplet_groups=len(results),
                   lt_peeled=res.n_peeled, lt_ge=res.n_ge, lt_missing=res.n_missing)
        if not res.ok:
            return DecodeResult(False, None, header, ext, {**rep, "error": "LT decode failure"}, cons_bodies)
        stored = Fountain.unsegment(res.segments, header.stored_len, header.global_seed)
    try:
        data = restore(stored, header.compressed)
    except Exception as exc:  # zstd raises its own error type
        return DecodeResult(False, None, header, ext, {**rep, "error": f"decompress: {exc}"}, cons_bodies)
    crc_ok = (zlib.crc32(data) & 0xFFFFFFFF) == header.crc32 and len(data) == header.original_len
    rep.update(crc32_ok=crc_ok, decode_s=time.perf_counter() - t0)
    return DecodeResult(crc_ok, data if crc_ok else None, header, ext, rep, cons_bodies)
