"""Noisy-read primitives: primer location/orientation, consensus, and clustering (edlib).

* :func:`locate_and_orient` finds F at the start and ``tail = revcomp(R)`` at the end of a
  read (or of its reverse complement) with edlib infix (``HW``) alignment, picks the
  orientation with the lower total primer edit distance, and returns the body between
  the primers. The body length varies with indels.
* :func:`consensus` aligns every read to a reference (edlib ``NW`` path), votes per
  reference column (base or deletion) and per gap (insertion), and iterates. The initial
  reference is the per-column majority of reads of exactly the expected length when there
  are enough of them, else a medoid. Column majority alone is *not* used as the answer:
  at a few % IDS, most exact-length reads still contain an insertion/deletion pair, whose
  shifted segment corrupts column votes. Measured: column-majority-only consensus at 3% IDS
  dropped from 0.965 exact at 10x to 0.735 at 20x.
* :func:`cluster` is a greedy edit-distance clustering, the fallback for reads whose index
  cannot be decoded.
* :func:`kmer_cluster` is the scalable grouping used on nanopore-like reads: greedy
  clustering where candidate centres come from shared k-mers and are verified with a
  bounded edlib distance. At ~10% read error a 15-nt index decodes exactly in only ~20% of
  reads (0.9^15), so grouping must not depend on the index field.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

import edlib
import numpy as np

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _find(query: str, text: str, k: int) -> tuple[int, int, int] | None:
    """Best infix match of ``query`` in ``text`` within ``k`` edits: (ed, start, end_inclusive)."""
    if not text:
        return None
    r = edlib.align(query, text, mode="HW", task="locations", k=k)
    if r["editDistance"] < 0 or not r["locations"]:
        return None
    start, end = r["locations"][0]
    return r["editDistance"], start, end


def locate_and_orient(read: str, forward: str, tail: str, max_primer_ed: int = 4,
                      margin: int = 80) -> tuple[str, bool, dict[str, Any]] | None:
    """Return ``(body, is_reverse, info)`` or ``None`` if the primers are not found.

    F is searched in the first ``len(F) + margin`` bases and ``tail`` in the last
    ``len(tail) + margin`` bases of each orientation. Both must match within
    ``max_primer_ed`` edits. Among the valid orientations, the one with the lower total edit
    distance wins (forward on ties).

    The default ``margin`` (80) covers sequencing-adapter/leader bases before F and after
    rc(R) in untrimmed nanopore reads (e.g. Dorado ``--no-trim``). It was 12, which
    rejected every flanked squigulator read (Phase 3 end-to-end run).
    """
    best = None
    for is_rev, s in ((False, read), (True, _rc(read))):
        hf = _find(forward, s[: len(forward) + margin], max_primer_ed)
        if hf is None:
            continue
        off = max(0, len(s) - len(tail) - margin)
        ht = _find(tail, s[off:], max_primer_ed)
        if ht is None:
            continue
        b0, b1 = hf[2] + 1, off + ht[1]
        if b1 < b0:
            continue
        tot = hf[0] + ht[0]
        if best is None or tot < best[0]:
            best = (tot, s[b0:b1], is_rev, {"ed_forward": hf[0], "ed_tail": ht[0]})
    if best is None:
        return None
    return best[1], best[2], best[3]


def _column_majority(bodies: list[str]) -> str:
    arr = np.frombuffer("".join(bodies).encode(), dtype=np.uint8).reshape(len(bodies), -1)
    out = []
    for col in arr.T:
        vals, cnt = np.unique(col, return_counts=True)
        out.append(chr(vals[np.argmax(cnt)]))
    return "".join(out)


def _medoid(bodies: list[str], rng: np.random.Generator, cap: int = 15) -> str:
    idx = np.arange(len(bodies))
    if len(idx) > cap:
        idx = rng.choice(idx, cap, replace=False)
    cand = [bodies[i] for i in idx]
    scores = [sum(edlib.align(a, b, task="distance")["editDistance"] for b in cand) for a in cand]
    return cand[int(np.argmin(scores))]


def _cigar_ops(cigar: str):
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            yield int(num), ch
            num = ""


def _align_vote(bodies: list[str], ref: str) -> str:
    """One round of alignment-to-reference voting."""
    L = len(ref)
    col_votes = [Counter() for _ in range(L)]  # base or "-" (deletion)
    ins_votes = [Counter() for _ in range(L + 1)]  # inserted string before column j ("" = none)
    for b in bodies:
        r = edlib.align(b, ref, mode="NW", task="path")
        qi = ti = 0
        ins_here: dict[int, str] = {}
        for n, op in _cigar_ops(r["cigar"]):
            if op in "=X":
                for _ in range(n):
                    col_votes[ti][b[qi]] += 1
                    qi += 1
                    ti += 1
            elif op == "I":  # extra bases in read, before ref column ti
                ins_here[ti] = ins_here.get(ti, "") + b[qi : qi + n]
                qi += n
            elif op == "D":  # ref bases missing from read
                for _ in range(n):
                    col_votes[ti]["-"] += 1
                    ti += 1
        for j in range(L + 1):
            ins_votes[j][ins_here.get(j, "")] += 1
    out = []
    for j in range(L + 1):
        ins, _ = ins_votes[j].most_common(1)[0]
        out.append(ins)
        if j < L:
            base, _ = col_votes[j].most_common(1)[0]
            if base != "-":
                out.append(base)
    return "".join(out)


def consensus(bodies: list[str], expected_len: int, rng: np.random.Generator | None = None,
              min_exact: int = 5, min_exact_frac: float = 0.3, rounds: int = 2) -> tuple[str, dict[str, Any]]:
    """Consensus body from a group of reads of the same oligo.

    Returns ``(sequence, info)`` with ``info = {"n", "n_exact", "method", "length"}``. The
    result is not padded or trimmed to ``expected_len``. Callers decide what to do with a
    wrong-length consensus (it normally fails CRC).
    """
    if not bodies:
        return "", {"n": 0, "n_exact": 0, "method": "none", "length": 0}
    rng = rng or np.random.default_rng(0)
    exact = [b for b in bodies if len(b) == expected_len]
    info: dict[str, Any] = {"n": len(bodies), "n_exact": len(exact)}
    if len(bodies) == 1:
        seq, info["method"] = bodies[0], "single"
    else:
        if len(exact) >= min_exact and len(exact) >= min_exact_frac * len(bodies):
            ref, info["method"] = _column_majority(exact), "align_vote(colmaj_init)"
        else:
            ref, info["method"] = _medoid(bodies, rng), "align_vote(medoid_init)"
        for _ in range(rounds):
            new = _align_vote(bodies, ref)
            if new == ref:
                break
            ref = new
        seq = ref
    info["length"] = len(seq)
    return seq, info


def cluster(bodies: list[str], max_ed_frac: float = 0.15, rng: np.random.Generator | None = None,
            max_centers: int | None = None) -> list[list[int]]:
    """Greedy clustering by edit distance (random visiting order).

    Each body joins the first existing centre within ``max_ed_frac * len(centre)`` edits
    (edlib with a ``k`` bound, so non-matches are cheap), otherwise it becomes a new
    centre. Returns lists of indices into ``bodies``.
    """
    rng = rng or np.random.default_rng(0)
    order = rng.permutation(len(bodies))
    centres: list[tuple[str, list[int]]] = []
    for i in order:
        b = bodies[i]
        placed = False
        for c, members in centres:
            k = max(1, int(max_ed_frac * len(c)))
            if abs(len(c) - len(b)) <= k and edlib.align(b, c, task="distance", k=k)["editDistance"] >= 0:
                members.append(int(i))
                placed = True
                break
        if not placed and (max_centers is None or len(centres) < max_centers):
            centres.append((b, [int(i)]))
    return [m for _, m in centres]


def kmer_cluster(bodies: list[str], k: int = 10, max_ed_frac: float = 0.3, min_votes: int = 2,
                 n_candidates: int = 3, order: list[int] | None = None) -> list[list[int]]:
    """Greedy clustering with k-mer candidate retrieval.

    Bodies are visited in ``order`` (default: as given). Each body is compared, with an edlib
    distance bounded by ``max_ed_frac * len``, against at most ``n_candidates`` existing
    centres that share at least ``min_votes`` distinct k-mers with it, and joins the closest
    one within the bound. Otherwise it becomes a new centre (a centre is its first member).

    ``max_ed_frac`` is a *read-to-read* bound. Measured (Phase 3, 160-nt bodies): reads of
    the same oligo after squigulator + Dorado hac differ by median 0.156, p99 0.287, max
    0.33 of the length; distinct oligo bodies differ by at least 0.36 (rotating ternary,
    median 0.425) and median 0.49 for whitened quaternary. 0.3 therefore splits ~1% of
    same-oligo pairs (harmless: duplicate groups resolve to the same index) rather than
    risk merging oligos. Merges are still survivable: the decoder also tries sub-groups by
    exact read index. Returns lists of indices into ``bodies``.
    """
    index: dict[str, list[int]] = {}
    centres: list[str] = []
    members: list[list[int]] = []
    for i in (order if order is not None else range(len(bodies))):
        b = bodies[i]
        kms = {b[j : j + k] for j in range(len(b) - k + 1)}
        votes: Counter = Counter()
        for km in kms:
            for c in index.get(km, ()):
                votes[c] += 1
        best, best_ed = -1, None
        for c, v in votes.most_common(n_candidates):
            if v < min_votes:
                break
            bound = int(max_ed_frac * max(len(b), len(centres[c])))
            ed = edlib.align(b, centres[c], task="distance", k=bound)["editDistance"]
            if ed >= 0 and (best_ed is None or ed < best_ed):
                best, best_ed = c, ed
        if best >= 0:
            members[best].append(int(i))
            continue
        cid = len(centres)
        centres.append(b)
        members.append([int(i)])
        for km in kms:
            index.setdefault(km, []).append(cid)
    return members
