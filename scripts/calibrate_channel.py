"""Channel calibration against real R10.4.1 reads of synthetic DNA-storage molecules.

Real data: Chen et al. 2025, Nat. Commun. "Approaching single-molecule assembly-free readout
from medium-length encoded DNA" (doi:10.1038/s41467-025-65004-7). Zenodo
10.5281/zenodo.16883332 (CC-BY-4.0): one POD5 (MinION FLO-MIN114, SQK-RBK114.24, 5 kHz) of
33 plasmids carrying LDPC/pseudo-noise-encoded DNA, plus the plasmid maps (SnapGene .dna).
Caveat: these are 6-43 kb plasmids read as ~2 kb transposase fragments, not a pool of
200-nt oligos. The encoded inserts are synthetic, but the plasmid backbone is natural vector
sequence. Both are reported (``region`` column).

Steps (``all`` runs them in order):
  real      map our own Dorado basecalls of the real POD5 (hac v6.0.0, sup v5.2.0)
  squig     simulate *the same reference intervals and strands* that the real hac reads cover
            with squigulator (dna-r10-min, 5 kHz) -> Dorado hac v6.0.0, then map
  badread   Badread v0.4.1 (nanopore2023 error + qscore model) on the plasmids, with identity
            matched to the real hac reads, no junk/random/chimera/adapter reads, then map
  report    per-base error by reference homopolymer run length, error-type mix, and the
            fraction of positions where >= 50% of covering reads err ("systematic sites", which a
            consensus keeps), at matched depth

Per-base error at reference base j: a substitution or deletion of j, or an insertion right
after j. Bases are grouped by the length of the reference run containing them (6 = 6+).
"""

from __future__ import annotations

import argparse
import gzip
import json
import struct
import subprocess
from collections import defaultdict
from pathlib import Path

import mappy
import numpy as np
import pandas as pd

from dnastore import nanopore
from dnastore.dnautil import revcomp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "real_r10" / "chen2025" / "Real-Time-Data-Readout-for-DNA-Storage"
POD5 = DATA / "real time readout data basecalling" / "pod5"
OUT = ROOT / "results" / "phase3" / "calibration"


def snapgene_seq(path: Path) -> str:
    b, i, seq = path.read_bytes(), 0, None
    while i < len(b):
        t, n = b[i], struct.unpack(">I", b[i + 1 : i + 5])[0]
        if t == 0:
            seq = b[i + 6 : i + 5 + n].decode().upper()
        i += 5 + n
    if seq is None:
        raise ValueError(f"no DNA packet in {path}")
    return seq


def load_refs() -> dict[str, str]:
    return {p.stem: snapgene_seq(p) for p in sorted((DATA / "Source data" / "Plasmid sequences").glob("*/*/*.dna"))}


def encoded_mask(refs: dict[str, str]) -> dict[str, np.ndarray]:
    """True where a reference base lies inside an exact copy of a published encoded sequence
    (Encoded sequences.txt), on either strand. The rest is vector/backbone."""
    enc = [l.strip().upper() for l in open(DATA / "Source data" / "Encoded sequences" / "Encoded sequences.txt")
           if l.strip()]
    probes = [e[i : i + 40] for e in enc for i in range(0, len(e) - 40, 40)]
    out = {}
    for name, s in refs.items():
        m = np.zeros(len(s), bool)
        for p in probes:
            for q in (p, revcomp(p)):
                j = s.find(q)
                while j >= 0:
                    m[j : j + 40] = True
                    j = s.find(q, j + 1)
        out[name] = m
    return out


def run_lengths(s: str) -> np.ndarray:
    out = np.zeros(len(s), np.int16)
    j = 0
    while j < len(s):
        k = j
        while k + 1 < len(s) and s[k + 1] == s[j]:
            k += 1
        out[j : k + 1] = k - j + 1
        j = k + 1
    return out


def read_fastq(path: Path) -> list[tuple[str, str]]:
    return [(rid, s) for rid, s, _ in nanopore.read_fastq(path)]


MIN_MAPQ = 20


def map_reads(reads: list[tuple[str, str]], refs: dict[str, str], min_mapq: int | None = None, min_len: int = 500):
    """Primary alignments: (read_id, ref, r_st, r_en, strand, per-ref-position error arrays)."""
    tmp = OUT / "refs.fa"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    tmp.write_text("".join(f">{k}\n{v}\n" for k, v in refs.items()))
    al = mappy.Aligner(str(tmp), preset="map-ont")
    min_mapq = MIN_MAPQ if min_mapq is None else min_mapq
    hits = []
    for rid, seq in reads:
        if len(seq) < min_len:
            continue
        for h in al.map(seq, cs=False):
            if not h.is_primary or h.mapq < min_mapq:
                continue
            ref = refs[h.ctg]
            q = seq if h.strand == 1 else revcomp(seq)
            qs = h.q_st if h.strand == 1 else len(seq) - h.q_en
            n = h.r_en - h.r_st
            sub, dele, ins = np.zeros(n, np.int8), np.zeros(n, np.int8), np.zeros(n, np.int8)
            ti, qi = 0, qs
            for length, op in h.cigar:
                if op == 0:  # M
                    a = np.frombuffer(q[qi : qi + length].encode(), np.uint8)
                    b = np.frombuffer(ref[h.r_st + ti : h.r_st + ti + length].encode(), np.uint8)
                    sub[ti : ti + length] = a != b
                    ti += length
                    qi += length
                elif op == 1:  # I (before ref position ti): count after base ti-1
                    if ti > 0:
                        ins[ti - 1] = 1
                    qi += length
                elif op == 2:  # D
                    dele[ti : ti + length] = 1
                    ti += length
            hits.append((rid, h.ctg, h.r_st, h.r_en, h.strand, sub, dele, ins))
            break
    return hits


def stats(hits, refs, masks, label: str, min_depth: int = 5) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rl = {k: run_lengths(v) for k, v in refs.items()}
    rows = []
    depth = {k: np.zeros(len(v), np.int32) for k, v in refs.items()}
    errc = {k: np.zeros(len(v), np.int32) for k, v in refs.items()}
    tot = defaultdict(int)
    for rid, ctg, a, b, strand, sub, dele, ins in hits:
        e = (sub | dele | ins).astype(np.int32)
        depth[ctg][a:b] += 1
        errc[ctg][a:b] += e
        tot["bases"] += b - a
        tot["sub"] += int(sub.sum())
        tot["del"] += int(dele.sum())
        tot["ins"] += int(ins.sum())
        r = np.minimum(rl[ctg][a:b], 6)
        reg = masks[ctg][a:b]
        for region, sel in (("encoded", reg), ("backbone", ~reg)):
            for L in range(1, 7):
                m = sel & (r == L)
                if m.any():
                    rows.append((label, region, L, int(m.sum()), int(e[m].sum())))
    df = pd.DataFrame(rows, columns=["channel", "region", "run_len", "bases", "errors"])
    per = df.groupby(["channel", "region", "run_len"], as_index=False)[["bases", "errors"]].sum()
    per["err_rate"] = per.errors / per.bases
    srows = []
    for ctg in refs:
        d, ec = depth[ctg], errc[ctg]
        ok = d >= min_depth
        r = np.minimum(rl[ctg], 6)
        for region, sel in (("encoded", masks[ctg]), ("backbone", ~masks[ctg])):
            for L in range(1, 7):
                m = ok & sel & (r == L)
                if m.any():
                    srows.append((label, region, L, int(m.sum()), int((ec[m] / d[m] >= 0.5).sum()),
                                  float(d[m].mean())))
    sdf = pd.DataFrame(srows, columns=["channel", "region", "run_len", "positions", "systematic", "mean_depth"])
    sdf = sdf.groupby(["channel", "region", "run_len"], as_index=False).agg(
        positions=("positions", "sum"), systematic=("systematic", "sum"), mean_depth=("mean_depth", "mean"))
    sdf["systematic_frac"] = sdf.systematic / sdf.positions
    b = max(1, tot["bases"])
    summ = {"channel": label, "mapped_reads": len(hits), "ref_bases": tot["bases"], "sub_rate": tot["sub"] / b,
            "del_rate": tot["del"] / b, "ins_rate": tot["ins"] / b,
            "err_rate": (tot["sub"] + tot["del"] + tot["ins"]) / b}
    return per, sdf, summ


def save_hits(hits, path: Path) -> None:
    pd.DataFrame([(h[0], h[1], h[2], h[3], h[4]) for h in hits],
                 columns=["read_id", "ref", "r_st", "r_en", "strand"]).to_csv(path, index=False)


def step_real(refs):
    out = {}
    for m in ("hac", "sup"):
        fq = OUT / "real_basecalls" / f"real_{m}.fastq"
        hits = map_reads(read_fastq(fq), refs)
        save_hits(hits, OUT / f"hits_real_{m}.csv")
        out[m] = hits
        print(f"real {m}: {len(hits)} mapped primary alignments", flush=True)
    return out


def step_squig(refs, profile: str = "dna-r10-min", seed: int = 11):
    iv = pd.read_csv(OUT / "hits_real_hac.csv")
    seqs = []
    for _, r in iv.iterrows():
        s = refs[r.ref][r.r_st : r.r_en]
        seqs.append(s if r.strand == 1 else revcomp(s))
    d = OUT / f"squig_{profile}"
    sim = nanopore.simulate_reads(seqs, [1] * len(seqs), d, seed=seed, profile=profile, rc_frac=0.0)
    pod5, _ = nanopore.to_pod5(sim["blow5"], d)
    fq, _ = nanopore.basecall(pod5, d, model="hac")
    for p in (sim["blow5"], pod5):
        Path(p).unlink(missing_ok=True)
    hits = map_reads(read_fastq(fq), refs)
    save_hits(hits, OUT / f"hits_squig_{profile}.csv")
    print(f"squigulator {profile}: {len(seqs)} intervals simulated, {len(hits)} mapped", flush=True)
    return hits


def step_badread(refs, identity: float, seed: int = 7):
    d = OUT / "badread"
    d.mkdir(parents=True, exist_ok=True)
    total = int(pd.read_csv(OUT / "hits_real_hac.csv").eval("r_en - r_st").sum())
    ref_len = sum(len(v) for v in refs.values())
    depth = max(1.0, total / ref_len)
    idp = f"{100 * identity:.2f},{min(99.9, 100 * identity + 4):.2f},{3.0}"
    cmd = [str(ROOT / ".venv" / "bin" / "badread"), "simulate", "--reference", str(OUT / "refs.fa"),
           "--quantity", f"{depth:.2f}x", "--error_model", "nanopore2023", "--qscore_model", "nanopore2023",
           "--identity", idp, "--length", "2000,1000", "--junk_reads", "0", "--random_reads", "0",
           "--chimeras", "0", "--glitches", "0,0,0", "--start_adapter_seq", "", "--end_adapter_seq", "",
           "--seed", str(seed)]
    with open(d / "badread.fastq", "w") as fh, open(d / "badread.log", "w") as lg:
        lg.write(" ".join(cmd) + "\n")
        subprocess.run(cmd, stdout=fh, stderr=lg, check=True)
    hits = map_reads(read_fastq(d / "badread.fastq"), refs)
    save_hits(hits, OUT / "hits_badread.csv")
    print(f"badread: identity {idp}, depth {depth:.2f}x, {len(hits)} mapped", flush=True)
    return hits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["all", "real", "squig", "badread", "report"], default="all", nargs="?")
    ap.add_argument("--min-mapq", type=int, default=20,
                    help="0 = sensitivity run: keeps reads tied between near-identical plasmids (MAPQ 0)")
    args = ap.parse_args()
    global MIN_MAPQ
    MIN_MAPQ = args.min_mapq
    tag = "" if args.min_mapq == 20 else f"_mapq{args.min_mapq}"
    OUT.mkdir(parents=True, exist_ok=True)
    refs = load_refs()
    masks = encoded_mask(refs)
    print(f"{len(refs)} plasmids, {sum(map(len, refs.values()))} nt, encoded fraction "
          f"{sum(m.sum() for m in masks.values()) / sum(map(len, refs.values())):.2f}", flush=True)
    chans = {}
    if args.step in ("all", "real", "report"):
        real = step_real(refs)
        chans["real_hac"], chans["real_sup"] = real["hac"], real["sup"]
    if args.step in ("all", "squig"):
        chans["squig_hac"] = step_squig(refs)
    if args.step in ("all", "badread"):
        _, _, s = stats(chans.get("real_hac") or step_real(refs)["hac"], refs, masks, "real_hac")
        chans["badread"] = step_badread(refs, 1 - s["err_rate"])
    if args.step == "report":
        for name, f in (("squig_hac", "hits_squig_dna-r10-min.csv"), ("badread", "hits_badread.csv")):
            if (OUT / f).exists():
                print(f"(report: re-mapping {name})", flush=True)
                fq = {"squig_hac": OUT / "squig_dna-r10-min" / "basecalls_hac.fastq",
                      "badread": OUT / "badread" / "badread.fastq"}[name]
                chans[name] = map_reads(read_fastq(fq), refs)
    pers, syss, summs = [], [], []
    for name, hits in chans.items():
        p, s, m = stats(hits, refs, masks, name)
        pers.append(p)
        syss.append(s)
        summs.append(m)
    per, sysd, summ = pd.concat(pers), pd.concat(syss), pd.DataFrame(summs)
    per.to_csv(OUT / f"err_by_runlen{tag}.csv", index=False)
    sysd.to_csv(OUT / f"systematic_by_runlen{tag}.csv", index=False)
    summ.to_csv(OUT / f"summary{tag}.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(summ.round(4).to_string(index=False))
        print(per.pivot_table(index=["region", "run_len"], columns="channel", values="err_rate").round(4))
        print(sysd.pivot_table(index=["region", "run_len"], columns="channel", values="systematic_frac").round(4))
        print(sysd.pivot_table(index=["region", "run_len"], columns="channel", values="positions"))
    (OUT / "run.json").write_text(json.dumps({"versions": nanopore.tool_versions(), "mappy": mappy.__version__,
                                              "dataset": "zenodo 10.5281/zenodo.16883332"}, indent=1, default=str))


if __name__ == "__main__":
    main()
