"""CRC-guided repair of a consensus body that failed its inner CRC.

Phase 3 found that squigulator + Dorado consensus bodies keep a sequence-determined error
floor (docs/notes/nanopore_pipeline.md), concentrated in homopolymer length errors. The repair
enumerates **homopolymer run-length edits of +-1** (a run of length 1 counts, so deleting or
duplicating a single base is included), at most ``max_edits`` = 2, such that the edited body
has the expected length. Each candidate is tested by the caller (index range, codec decode,
CRC-32 bound to the index), and the first one that passes is accepted.

**Order and budget.** The candidates are ranked by read support. For every run of the
consensus, the fraction of aligned reads that show the run one shorter (any deletion inside
it) or one longer (an insertion of the run's base inside or next to it) is counted, plus
a pseudo-count, so unsupported edits are still tried, but last. A pair's score is the
product of its two single-edit scores. At most ``budget`` candidates are tested. The CRC
false-accept probability is then at most ``budget * 2**-32`` per repaired oligo.

Substitutions are *not* enumerated, so the residual error after repair is reported.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterator

import edlib

from .consensus import _cigar_ops


def runs(seq: str) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for ch in seq:
        if out and out[-1][0] == ch:
            out[-1] = (ch, out[-1][1] + 1)
        else:
            out.append((ch, 1))
    return out


def run_support(bodies: list[str], ref: str) -> list[tuple[float, float]]:
    """Per run of ``ref``: (fraction of reads supporting length - 1, ... length + 1)."""
    rr = runs(ref)
    starts, s = [], 0
    for _, n in rr:
        starts.append(s)
        s += n
    run_of = [i for i, (_, n) in enumerate(rr) for _ in range(n)]
    minus = [0] * len(rr)
    plus = [0] * len(rr)
    n_reads = 0
    for b in bodies:
        if not b:
            continue
        n_reads += 1
        r = edlib.align(b, ref, mode="NW", task="path")
        qi = ti = 0
        hit_m: set[int] = set()
        hit_p: set[int] = set()
        for n, op in _cigar_ops(r["cigar"]):
            if op in "=X":
                qi += n
                ti += n
            elif op == "D":
                for t in range(ti, ti + n):
                    hit_m.add(run_of[t])
                ti += n
            elif op == "I":
                ins = b[qi : qi + n]
                # the insertion sits before ref column ti: next to run_of[ti-1] and run_of[ti]
                for t in (ti - 1, ti):
                    if 0 <= t < len(ref) and set(ins) == {ref[t]}:
                        hit_p.add(run_of[t])
                qi += n
        for i in hit_m:
            minus[i] += 1
        for i in hit_p:
            plus[i] += 1
    if n_reads == 0:
        return [(0.0, 0.0)] * len(rr)
    return [(m / n_reads, p / n_reads) for m, p in zip(minus, plus)]


def _apply(rr: list[tuple[str, int]], edits: dict[int, int]) -> str:
    return "".join(ch * (n + edits.get(i, 0)) for i, (ch, n) in enumerate(rr))


def candidates(ref: str, target_len: int, support: list[tuple[float, float]] | None = None,
               max_edits: int = 2, budget: int = 2000, pseudo: float = 0.02) -> Iterator[tuple[str, int]]:
    """Yield ``(candidate, n_edits)`` with ``len == target_len``, best-supported first."""
    rr = runs(ref)
    delta = target_len - len(ref)
    if abs(delta) > max_edits or delta == 0 and max_edits < 2:
        return
    sup = support or [(0.0, 0.0)] * len(rr)
    singles = []  # (score, run, +-1)
    for i, (_, n) in enumerate(rr):
        singles.append((sup[i][0] + pseudo, i, -1))
        singles.append((sup[i][1] + pseudo, i, +1))
    singles.sort(reverse=True)
    scored: list[tuple[float, dict[int, int]]] = []
    if abs(delta) == 1:
        scored = [(s, {i: d}) for s, i, d in singles if d == delta]
    elif max_edits >= 2:
        want = (-1, +1) if delta == 0 else ((delta // 2, delta // 2))
        a = [x for x in singles if x[2] == want[0]]
        b = [x for x in singles if x[2] == want[1]]
        for sa, ia, da in a:
            for sb, ib, db in b:
                if delta == 0 and ia == ib:
                    continue
                if delta != 0 and ib < ia:
                    continue  # unordered pairs of same-sign edits
                e: dict[int, int] = {ia: da}
                e[ib] = e.get(ib, 0) + db
                if any(rr[i][1] + v < 0 for i, v in e.items()):
                    continue
                scored.append((sa * sb, e))
        scored.sort(key=lambda x: -x[0])
    seen: set[str] = set()
    for _, e in scored:
        if len(seen) >= budget:
            break
        c = _apply(rr, e)
        if c in seen or len(c) != target_len:
            continue
        seen.add(c)
        yield c, sum(abs(v) for v in e.values())


@dataclass
class RepairResult:
    body: str | None
    tried: int
    edits: int


def repair(ref: str, bodies: list[str], target_len: int, accept: Callable[[str], bool], max_edits: int = 2,
           budget: int = 2000) -> RepairResult:
    """Return the first candidate for which ``accept(candidate)`` is true."""
    if abs(target_len - len(ref)) > max_edits:
        return RepairResult(None, 0, 0)
    sup = run_support(bodies, ref) if bodies else None
    tried = 0
    for cand, ne in candidates(ref, target_len, sup, max_edits, budget):
        tried += 1
        if accept(cand):
            return RepairResult(cand, tried, ne)
    return RepairResult(None, tried, 0)
