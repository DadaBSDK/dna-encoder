// DNASTORE-RS-1: error-tolerant DNA archive for the browser.
//
// The file is cut into 32-byte payloads. Each strand is
//     [index: 4 bytes][payload: 32 bytes][CRC-32(index || payload): 4 bytes]   -> 160 nt
// The index is masked and payload+CRC are XOR-whitened with a keystream keyed by the index,
// so runs of identical bytes do not become long homopolymers. Strands are written as FASTA
// records, but the decoder ignores record names and order: identity comes from the DNA.
//
// Redundancy and recovery, in order:
//   1. CRC-32 rejects any damaged strand, so damage stays local to one strand.
//   2. A rejected strand with one substitution, insertion or deletion is repaired by trying
//      every single edit and keeping the one that passes the CRC.
//   3. Systematic Reed-Solomon over GF(256) across strands (k data + p parity per block,
//      p = ceil(k * parity_permille / 1000)) rebuilds any p missing strands per block.
//   4. Compromise: if a block has more missing strands than parity, its data is still
//      returned. Damaged strands contribute their best guess and absent ones zeros; those
//      bytes are reported as uncertain. Bytes NOT reported uncertain are exact.
// The metadata (name, size, SHA-256, parity) is stored in 8 copies of dedicated header strands.

import { Sha256 } from './stream.js';

export const FORMAT = 'DNASTORE-RS-1';
export const PAYLOAD = 32;
export const STRAND_BYTES = 4 + PAYLOAD + 4;
export const STRAND_NT = STRAND_BYTES * 4;
export const HEADER_BASE = 0xffffff00;
export const HEADER_COPIES = 8;
export const MAX_BYTES = 64 * 1024 * 1024;
const INDEX_MASK = 0x9e3779b9;

export class ResilientError extends Error {}

// --------------------------------------------------------------------------- CRC-32, whitening

const CRC_TABLE = new Uint32Array(256).map((_, n) => {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  return c >>> 0;
});

export function crc32(bytes, start = 0, end = bytes.length) {
  let c = 0xffffffff;
  for (let i = start; i < end; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

function whiten(buf, index) {
  let x = (Math.imul(index ^ 0x5bd1e995, 0x01000193) ^ 0x2545f491) >>> 0 || 1;
  for (let i = 4; i < STRAND_BYTES; i++) {
    x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0;
    buf[i] ^= x & 0xff;
  }
}

// --------------------------------------------------------------------------- strand <-> DNA

const BASES = 'ACGT';
const BASE_VALUE = { A: 0, C: 1, G: 2, T: 3 };

function strandBytes(index, payload) {
  const buf = new Uint8Array(STRAND_BYTES);
  new DataView(buf.buffer).setUint32(0, index >>> 0);
  buf.set(payload, 4);
  new DataView(buf.buffer).setUint32(4 + PAYLOAD, crc32(buf, 0, 4 + PAYLOAD));
  whiten(buf, index);
  new DataView(buf.buffer).setUint32(0, (index ^ INDEX_MASK) >>> 0);
  return buf;
}

function bytesToDna(buf) {
  let s = '';
  for (const b of buf) s += BASES[b >> 6] + BASES[(b >> 4) & 3] + BASES[(b >> 2) & 3] + BASES[b & 3];
  return s;
}

const work = new Uint8Array(STRAND_BYTES);

/** Check packed strand bytes ``raw`` (40 bytes as written). Returns {index, payload, ok}. */
function checkRaw(raw) {
  work.set(raw);
  const index = (((work[0] << 24) | (work[1] << 16) | (work[2] << 8) | work[3]) ^ INDEX_MASK) >>> 0;
  work[0] = index >>> 24; work[1] = index >>> 16; work[2] = index >>> 8; work[3] = index;
  whiten(work, index);
  const tag = ((work[36] << 24) | (work[37] << 16) | (work[38] << 8) | work[39]) >>> 0;
  return { index, ok: tag === crc32(work, 0, 4 + PAYLOAD), payload: work.slice(4, 4 + PAYLOAD) };
}

/** Base values (0..3, or 4 for anything else) of a sequence. */
const CODE_VALUE = new Uint8Array(128).fill(4);
CODE_VALUE[65] = 0; CODE_VALUE[67] = 1; CODE_VALUE[71] = 2; CODE_VALUE[84] = 3;
function baseValues(seq) {
  const v = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) { const c = seq.charCodeAt(i); v[i] = c < 128 ? CODE_VALUE[c] : 4; }
  return v;
}

/** Pack exactly STRAND_NT base values into bytes; null if any base is not A/C/G/T. */
function pack(vals) {
  const raw = new Uint8Array(STRAND_BYTES);
  for (let i = 0; i < STRAND_NT; i++) {
    if (vals[i] > 3) return null;
    raw[i >> 2] |= vals[i] << (6 - 2 * (i & 3));
  }
  return raw;
}

/** Parse one strand. ``ok`` says whether the CRC passed; null if it cannot be read at all. */
function parseStrand(seq) {
  if (seq.length !== STRAND_NT) return null;
  const raw = pack(baseValues(seq));
  return raw && checkRaw(raw);
}

// CRC-32 is affine: for equal-length messages crc(m ^ e) = crc(m) ^ crc(e) ^ crc(0).
// DELTA[j * 256 + x] = crc(e) ^ crc(0) where e is x at byte j of a 36-byte message, so a
// one-byte change is tested with one table lookup instead of a full CRC.
const CRC_LEN = 4 + PAYLOAD;
const DELTA = (() => {
  const t = new Uint32Array(CRC_LEN * 256), m = new Uint8Array(CRC_LEN), zero = crc32(m);
  for (let j = 0; j < CRC_LEN; j++) for (let x = 1; x < 256; x++) { m[j] = x; t[j * 256 + x] = crc32(m) ^ zero; m[j] = 0; }
  return t;
})();

function packAny(vals) {
  const raw = new Uint8Array(STRAND_BYTES);
  for (let i = 0; i < STRAND_NT; i++) raw[i >> 2] |= (vals[i] & 3) << (6 - 2 * (i & 3));
  return raw;
}

/** Try every single-base substitution of an exact-length read (positions limited to ``only``
 * when it is set, e.g. the position of an N). Payload and tag positions cost one lookup each. */
function repairSubstitution(vals, only = null, withIndex = true) {
  const raw = packAny(vals);
  const positions = only ?? Array.from({ length: STRAND_NT }, (_, i) => i);
  // index bytes change the whitening keystream: full check
  for (const i of positions) {
    if (i >= 16 || !withIndex) continue;
    const j = i >> 2, sh = 6 - 2 * (i & 3), keep = raw[j];
    for (let b = 0; b < 4; b++) {
      if (b === vals[i]) continue;
      raw[j] = (keep & ~(3 << sh)) | (b << sh);
      const r = checkRaw(raw);
      if (r.ok) return r;
    }
    raw[j] = keep;
  }
  // payload/tag bytes: dewhiten once, then CRC deltas
  const base = checkRaw(raw);
  const plain = work.slice();                      // dewhitened bytes from checkRaw
  const c0 = crc32(plain, 0, CRC_LEN);
  const tag = ((plain[36] << 24) | (plain[37] << 16) | (plain[38] << 8) | plain[39]) >>> 0;
  for (const i of positions) {
    if (i < 16) continue;
    const j = i >> 2, sh = 6 - 2 * (i & 3), cur = (raw[j] >> sh) & 3;
    for (let b = 0; b < 4; b++) {
      if (b === cur) continue;
      const x = ((b ^ cur) << sh) & 0xff;
      const hit = j < CRC_LEN ? (c0 ^ DELTA[j * 256 + x]) >>> 0 === tag
        : c0 === (tag ^ (x << (8 * (39 - j)))) >>> 0;
      if (hit) {
        raw[j] ^= x;
        const r = checkRaw(raw);
        if (r.ok) return r;
        raw[j] ^= x;
      }
    }
  }
  return base.ok ? base : null;
}

/** CRC-guided repair of one substitution, insertion or deletion (or one N). */
/** ``withIndex = false`` skips edits inside the 16-nt index field, which need full checks;
 * the decoder only pays for them when parity cannot cover the remaining gaps. */
function repairStrand(vals, withIndex = true) {
  const n = vals.length;
  if (n === STRAND_NT) {
    const bad = [];
    for (let i = 0; i < n; i++) if (vals[i] > 3) bad.push(i);
    if (bad.length > 1) return null;
    return repairSubstitution(vals, bad.length ? bad : null, withIndex);
  }
  if (n !== STRAND_NT - 1 && n !== STRAND_NT + 1) return null;
  for (let i = 0; i < n; i++) if (vals[i] > 3) return null;
  const cand = new Uint8Array(STRAND_NT);
  // candidate with the edit at position i: deletion repair inserts base b at i (read is
  // one short); insertion repair drops read[i] (read is one long)
  const build = (i, b) => {
    if (n < STRAND_NT) { cand.set(vals.subarray(0, i)); cand[i] = b; cand.set(vals.subarray(i), i + 1); }
    else { cand.set(vals.subarray(0, i)); cand.set(vals.subarray(i + 1), i); }
  };
  const tries = n < STRAND_NT ? [0, 1, 2, 3] : [-1];
  for (const b of tries) {
    // edits inside the index field change the whitening: full checks
    for (let i = 0; i < 16 && withIndex; i++) {
      build(i, b);
      const r = checkRaw(packAny(cand));
      if (r.ok) return r;
    }
    // from i = 16 on, the index is fixed and consecutive candidates differ in one or two
    // bases: walk them, updating the dewhitened CRC and tag with table lookups
    build(16, b);
    const raw = packAny(cand);
    checkRaw(raw);
    let c = crc32(work, 0, CRC_LEN);
    let t = ((work[36] << 24) | (work[37] << 16) | (work[38] << 8) | work[39]) >>> 0;
    const setBase = (pos, v) => {
      const j = pos >> 2, sh = 6 - 2 * (pos & 3), x = ((((raw[j] >> sh) & 3) ^ v) << sh) & 0xff;
      if (!x) return;
      raw[j] ^= x;
      if (j < CRC_LEN) c = (c ^ DELTA[j * 256 + x]) >>> 0; else t = (t ^ (x << (8 * (39 - j)))) >>> 0;
    };
    const last = n < STRAND_NT ? n : n - 1;
    for (let i = 16; ; i++) {
      if (c === t) { const r = checkRaw(raw); if (r.ok) return r; }
      if (i >= last) break;
      if (n < STRAND_NT) { setBase(i, vals[i]); setBase(i + 1, b); } else setBase(i, vals[i]);
    }
  }
  return null;
}

/** True if ``vals`` is within a few edits of the verified strand (index, payload).
 * Whitening makes different strands differ in ~75% of positions, copies in a few. */
function nearCopy(vals, payload, index) {
  const ref = baseValues(bytesToDna(strandBytes(index, payload)));
  const n = Math.min(ref.length, vals.length), half = n >> 1;
  let diff = Math.abs(ref.length - vals.length);
  for (let i = 0; i < half; i++) diff += ref[i] !== vals[i];
  for (let i = 1; i <= n - half; i++) diff += ref[ref.length - i] !== vals[vals.length - i];
  return diff <= 40;
}

/** Best-guess payload of a read that cannot be verified: one with an indel is cut or
 * padded at its end, so bytes before the indel are still right. Only used for bytes that
 * are reported as uncertain. */
function guessOf(vals) {
  if (Math.abs(vals.length - STRAND_NT) > 8 || vals.length < 16) return null;
  const fixed = new Uint8Array(STRAND_NT);
  fixed.set(vals.subarray(0, STRAND_NT));
  const r = checkRaw(packAny(fixed));
  return r.payload;
}

/** Index a read claims (from its first 16 bases), without verification; -1 if unreadable. */
function claimedIndex(vals) {
  if (vals.length < 16) return -1;
  let v = 0;
  for (let i = 0; i < 16; i++) { if (vals[i] > 3) return -1; v = (v * 4) + vals[i]; }
  return (v ^ INDEX_MASK) >>> 0;
}

const BAND = 10;
const INF = 0xffff;
let dpBuf = new Uint16Array(0);

/** Banded global alignment of ``read`` to ``ref`` (unit costs, |i - j| <= BAND). Returns
 * per-ref-position votes: base value (0..3), 5 = deleted; and inserted bases before each
 * ref position. Reads of one strand differ by a few indels, so the band is ample. */
function alignVotes(ref, read) {
  const n = ref.length, m = read.length, w = m + 1;
  if ((n + 1) * w > dpBuf.length) dpBuf = new Uint16Array((n + 1) * w);
  const d = dpBuf;
  for (let i = 0; i <= n; i++) {
    const lo = Math.max(0, i - BAND - 1), hi = Math.min(m, i + BAND + 1);
    for (let j = lo; j <= hi; j++) {
      if (i === 0) { d[j] = j <= BAND ? j : INF; continue; }
      if (j === 0) { d[i * w] = i <= BAND ? i : INF; continue; }
      if (Math.abs(i - j) > BAND) { d[i * w + j] = INF; continue; }
      const sub = d[(i - 1) * w + j - 1] + (ref[i - 1] === read[j - 1] ? 0 : 1);
      const up = d[(i - 1) * w + j] + 1, left = d[i * w + j - 1] + 1;
      d[i * w + j] = Math.min(sub, up, left, INF);
    }
  }
  const col = new Uint8Array(n).fill(5), insBefore = Array.from({ length: n + 1 }, () => []);
  if (d[n * w + m] >= INF || Math.abs(n - m) > BAND) return { col, insBefore };
  let i = n, j = m;
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && d[i * w + j] === d[(i - 1) * w + j - 1] + (ref[i - 1] === read[j - 1] ? 0 : 1)) { col[i - 1] = read[j - 1]; i--; j--; }
    else if (i > 0 && d[i * w + j] === d[(i - 1) * w + j] + 1) { col[i - 1] = 5; i--; }
    else { insBefore[i].push(read[j - 1]); j--; }
  }
  return { col, insBefore };
}

/** Majority consensus of several reads of one strand, refined twice. */
function consensus(reads) {
  let ref = reads.find((r) => r.length === STRAND_NT) || reads[0];
  for (let round = 0; round < 2; round++) {
    const n = ref.length;
    const votes = Array.from({ length: n }, () => new Uint16Array(6));
    const ins = Array.from({ length: n + 1 }, () => new Uint16Array(4));
    const insCount = new Uint16Array(n + 1);
    for (const r of reads) {
      const { col, insBefore } = alignVotes(ref, r);
      for (let i = 0; i < n; i++) votes[i][col[i]]++;
      for (let i = 0; i <= n; i++) if (insBefore[i].length) { insCount[i]++; ins[i][insBefore[i][insBefore[i].length - 1] & 3]++; }
    }
    const out = [];
    for (let i = 0; i <= n; i++) {
      if (insCount[i] * 2 > reads.length) out.push(ins[i].indexOf(Math.max(...ins[i])));
      if (i === n) break;
      const v = votes[i];
      let best = 0;
      for (let b = 1; b < 4; b++) if (v[b] > v[best]) best = b;
      if (v[5] <= v[best]) out.push(best);
    }
    ref = Uint8Array.from(out);
  }
  return ref;
}

// --------------------------------------------------------------------------- GF(256) Reed-Solomon

const EXP = new Uint8Array(512), LOG = new Uint8Array(256);
for (let i = 0, x = 1; i < 255; i++) {
  EXP[i] = x; LOG[x] = i;
  x <<= 1; if (x & 0x100) x ^= 0x11d;
}
for (let i = 255; i < 512; i++) EXP[i] = EXP[i - 255];
const mul = (a, b) => (a && b ? EXP[LOG[a] + LOG[b]] : 0);
const inv = (a) => EXP[255 - LOG[a]];

/** Inverse of a square matrix over GF(256) (array of Uint8Array rows). */
function invert(m) {
  const n = m.length;
  const a = m.map((row, i) => { const r = new Uint8Array(2 * n); r.set(row); r[n + i] = 1; return r; });
  for (let c = 0; c < n; c++) {
    let p = c;
    while (p < n && !a[p][c]) p++;
    if (p === n) throw new ResilientError('singular RS matrix');
    [a[c], a[p]] = [a[p], a[c]];
    const f = inv(a[c][c]);
    for (let j = 0; j < 2 * n; j++) a[c][j] = mul(a[c][j], f);
    for (let r = 0; r < n; r++) {
      if (r === c || !a[r][c]) continue;
      const g = a[r][c];
      for (let j = 0; j < 2 * n; j++) a[r][j] ^= mul(g, a[c][j]);
    }
  }
  return a.map((row) => row.slice(n));
}

/** Fill the rows listed in ``erased`` so that every column is an RS codeword.
 * ``rows``: n Uint8Array(PAYLOAD) (erased rows' contents ignored); parity check
 * H[i][j] = alpha^(i * (n - 1 - j)), i < p. Requires erased.length <= p. */
function solveErasures(rows, erased, p) {
  const n = rows.length, e = erased.length;
  if (!e) return;
  const isErased = new Uint8Array(n);
  for (const j of erased) isErased[j] = 1;
  const h = (i, j) => EXP[(i * (n - 1 - j)) % 255];
  const minv = invert(Array.from({ length: e }, (_, i) => Uint8Array.from(erased, (j) => h(i, j))));
  const synd = Array.from({ length: e }, () => new Uint8Array(PAYLOAD));
  for (let i = 0; i < e; i++) {
    for (let j = 0; j < n; j++) {
      if (isErased[j]) continue;
      const coef = h(i, j), row = rows[j], s = synd[i];
      if (!coef) continue;
      for (let c = 0; c < PAYLOAD; c++) s[c] ^= mul(coef, row[c]);
    }
  }
  erased.forEach((j, r) => {
    const out = new Uint8Array(PAYLOAD);
    for (let i = 0; i < e; i++) {
      const coef = minv[r][i];
      if (!coef) continue;
      for (let c = 0; c < PAYLOAD; c++) out[c] ^= mul(coef, synd[i][c]);
    }
    rows[j] = out;
  });
}

// --------------------------------------------------------------------------- layout

export function parityRows(k, permille) { return Math.ceil((k * permille) / 1000); }

/** Balanced RS blocks [[k, p], ...] for ``nData`` data strands. */
export function planBlocks(nData, permille) {
  if (!nData) return [];
  let kmax = 255;
  while (kmax > 1 && kmax + parityRows(kmax, permille) > 255) kmax--;
  const nb = Math.ceil(nData / kmax);
  const base = Math.floor(nData / nb), extra = nData % nb;
  return Array.from({ length: nb }, (_, i) => { const k = base + (i < extra ? 1 : 0); return [k, parityRows(k, permille)]; });
}

// --------------------------------------------------------------------------- encode

/** Encode bytes into a FASTA archive string. ``permille`` = parity strands per 1000 data strands. */
export async function encodeResilient(blob, name, permille = 250) {
  if (!Number.isInteger(permille) || permille < 0 || permille > 4000) throw new ResilientError('parity must be 0..4000 permille');
  if (blob.size > MAX_BYTES) throw new ResilientError(`error-tolerant mode is limited to ${MAX_BYTES / 1048576} MB in the browser`);
  const data = new Uint8Array(await blob.arrayBuffer());
  const sha256 = new Sha256().update(data).hex();
  const meta = JSON.stringify({ format: FORMAT, name: String(name), bytes: data.length, sha256, parity: permille });
  const metaBytes = new TextEncoder().encode(meta);
  const hdr = new Uint8Array(2 + metaBytes.length);
  hdr[0] = metaBytes.length >> 8; hdr[1] = metaBytes.length & 0xff; hdr.set(metaBytes, 2);
  const nFrag = Math.ceil(hdr.length / PAYLOAD);
  if (metaBytes.length > 0xffff || nFrag > 255) throw new ResilientError('filename too long');

  const out = [];
  const add = (index, payload) => { out.push(`>s${out.length}\n${bytesToDna(strandBytes(index, payload))}\n`); };
  for (let copy = 0; copy < HEADER_COPIES; copy++) {
    for (let f = 0; f < nFrag; f++) {
      const p = new Uint8Array(PAYLOAD);
      p.set(hdr.subarray(f * PAYLOAD, (f + 1) * PAYLOAD));
      add(HEADER_BASE + f, p);
    }
  }
  const blocks = planBlocks(Math.ceil(data.length / PAYLOAD), permille);
  let q = 0, d = 0;
  for (const [k, p] of blocks) {
    const rows = [];
    for (let i = 0; i < k; i++, d++) {
      const r = new Uint8Array(PAYLOAD);
      r.set(data.subarray(d * PAYLOAD, (d + 1) * PAYLOAD));
      rows.push(r);
    }
    for (let i = 0; i < p; i++) rows.push(new Uint8Array(PAYLOAD));
    solveErasures(rows, Array.from({ length: p }, (_, i) => k + i), p);
    for (const r of rows) add(q++, r);
  }
  const text = out.join('');
  return {
    archive: new Blob([text], { type: 'text/plain' }), text,
    result: { name: String(name), bytes: data.length, sha256, parity: permille, strands: out.length,
      dataStrands: blocks.reduce((s, b) => s + b[0], 0), parityStrands: blocks.reduce((s, b) => s + b[1], 0),
      headerStrands: nFrag * HEADER_COPIES, nucleotides: out.length * STRAND_NT },
  };
}

// --------------------------------------------------------------------------- decode

/** Candidate strand sequences from FASTA or one-sequence-per-line text. Names are ignored. */
function candidates(text) {
  const seqs = [];
  const fasta = /^\s*>/.test(text);
  let cur = '';
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line[0] === '>') { if (cur) seqs.push(cur); cur = ''; continue; }
    if (line[0] === ';') continue;
    cur += line.toUpperCase().replace(/\s+/g, '');
    if (!fasta) { seqs.push(cur); cur = ''; }
  }
  if (cur) seqs.push(cur);
  return seqs;
}

export function looksResilient(text) {
  return /^\s*>/.test(text) || /^[ACGTNacgtn\s]+$/.test(text.slice(0, 2000));
}

/** Metadata from the verified header strands; null (or a ResilientError if ``strict``). */
function readMeta(valid, strict = false) {
  const fail = (msg) => { if (strict) throw new ResilientError(msg); return null; };
  const f0 = valid.get(HEADER_BASE);
  if (!f0) return fail('metadata strands lost: every copy of the header was damaged beyond repair');
  const len = (f0[0] << 8) | f0[1];
  const nFrag = Math.ceil((2 + len) / PAYLOAD);
  const hdr = new Uint8Array(nFrag * PAYLOAD);
  for (let f = 0; f < nFrag; f++) {
    const frag = valid.get(HEADER_BASE + f);
    if (!frag) return fail('metadata strands lost: a header fragment was damaged in every copy');
    hdr.set(frag, f * PAYLOAD);
  }
  let meta;
  try { meta = JSON.parse(new TextDecoder('utf-8', { fatal: true }).decode(hdr.subarray(2, 2 + len))); } catch { meta = null; }
  if (!meta || meta.format !== FORMAT || !Number.isSafeInteger(meta.bytes) || meta.bytes < 0 || typeof meta.sha256 !== 'string'
      || !Number.isInteger(meta.parity) || meta.parity < 0 || meta.parity > 4000) {
    return fail('not a DNASTORE-RS-1 archive (metadata unreadable)');
  }
  if (meta.bytes > MAX_BYTES) return fail('archive declares a file larger than this page supports');
  return meta;
}

/** Decode with redundancy and compromise. Never throws for damaged DNA: returns a report.
 * Throws only if no metadata (header) strand survives, since the file size is then unknown. */
export function decodeResilient(text) {
  const report = { candidates: 0, intact: 0, repaired: 0, consensus: 0, rejected: 0, duplicates: 0, skipped: 0 };
  const valid = new Map();       // index -> payload (CRC-verified)
  const guesses = new Map();     // index -> payload (failed CRC, best guess)
  const accept = (r) => { if (valid.has(r.index)) report.duplicates++; else valid.set(r.index, r.payload); };

  // pass 1: reads that verify as they are
  const pending = [];
  for (const seq of candidates(text)) {
    report.candidates++;
    const vals = baseValues(seq);
    const raw = vals.length === STRAND_NT ? pack(vals) : null;
    const r = raw && checkRaw(raw);
    if (r && r.ok) { report.intact++; accept(r); continue; }
    pending.push({ vals, claim: r ? r.index : claimedIndex(vals), guess: r ? r.payload : guessOf(vals) });
  }
  // pass 2: single-edit repair outside the index field. A read whose claimed index is
  // already verified is skipped only if it really is a copy of that strand: an error in the
  // index field makes a read claim a neighbour's index.
  let failed = [];
  const isCopy = (p) => valid.has(p.claim) && nearCopy(p.vals, valid.get(p.claim), p.claim);
  for (const p of pending) {
    if (isCopy(p)) { report.skipped++; continue; }
    const r = repairStrand(p.vals, false);
    if (r) { report.repaired++; accept(r); } else failed.push(p);
  }

  // Passes 3-4 are expensive; run them only if the header or some RS block is still short.
  const short = () => {
    const meta = readMeta(valid);
    if (!meta) return true;
    let q = 0;
    for (const [k, p] of planBlocks(Math.ceil(meta.bytes / PAYLOAD), meta.parity)) {
      let e = 0;
      for (let j = 0; j < k + p; j++) e += !valid.has(q + j);
      if (e > p) return true;
      q += k + p;
    }
    return false;
  };
  if (short()) {
    // pass 3: edits inside the index field
    const still = [];
    for (const p of failed) {
      if (isCopy(p)) { report.skipped++; continue; }
      const r = repairStrand(p.vals, true);
      if (r) { report.repaired++; accept(r); } else still.push(p);
    }
    failed = still;
    // pass 4: consensus of the still-failing reads that claim the same missing index
    const groups = new Map();
    for (const p of failed) {
      if (isCopy(p)) { report.skipped++; continue; }
      if (!groups.has(p.claim)) groups.set(p.claim, []);
      groups.get(p.claim).push(p);
    }
    for (const [claim, group] of groups) {
      if (group.length >= 2) {
        const cons = consensus(group.map((p) => p.vals));
        const r = cons.length === STRAND_NT ? checkRaw(packAny(cons)) : null;
        const ok = r && r.ok ? r : repairStrand(cons, true);
        if (ok && !valid.has(ok.index)) { report.consensus++; accept(ok); continue; }
        if (r && r.index === claim && !guesses.has(claim)) guesses.set(claim, r.payload);
      }
      report.rejected += group.length;
      for (const p of group) if (p.guess && !guesses.has(claim)) guesses.set(claim, p.guess);
    }
  } else {
    report.rejected += failed.filter((p) => !isCopy(p)).length;
  }

  const meta = readMeta(valid, true);

  // data
  const blocks = planBlocks(Math.ceil(meta.bytes / PAYLOAD), meta.parity);
  const out = new Uint8Array(blocks.reduce((s, b) => s + b[0], 0) * PAYLOAD);
  const uncertain = [];          // [start, end) byte ranges in the output
  let q = 0, d = 0, rebuilt = 0, guessed = 0, missing = 0, failedBlocks = 0;
  for (const [k, p] of blocks) {
    const n = k + p;
    const rows = [], erased = [];
    for (let j = 0; j < n; j++) {
      const v = valid.get(q + j);
      rows.push(v ? v.slice() : new Uint8Array(PAYLOAD));
      if (!v) erased.push(j);
    }
    if (erased.length <= p) {
      rebuilt += erased.filter((j) => j < k).length;
      solveErasures(rows, erased, p);
    } else {
      failedBlocks++;
      for (const j of erased) {
        if (j >= k) continue;
        const g = guesses.get(q + j);
        if (g) { rows[j] = g.slice(); guessed++; } else missing++;
        uncertain.push([(d + j) * PAYLOAD, (d + j + 1) * PAYLOAD]);
      }
    }
    for (let j = 0; j < k; j++) out.set(rows[j], (d + j) * PAYLOAD);
    q += n; d += k;
  }
  const data = out.subarray(0, meta.bytes);
  const ranges = [];
  for (const [s, e] of uncertain) {
    const a = Math.min(s, meta.bytes), b = Math.min(e, meta.bytes);
    if (a >= b) continue;
    if (ranges.length && ranges[ranges.length - 1][1] === a) ranges[ranges.length - 1][1] = b; else ranges.push([a, b]);
  }
  const uncertainBytes = ranges.reduce((s, [a, b]) => s + b - a, 0);
  const sha = new Sha256().update(data).hex();
  return {
    file: new Blob([data]), bytes: data,
    verified: sha === meta.sha256, meta,
    report: { ...report, dataStrands: d, rebuiltByParity: rebuilt, guessed, missing, failedBlocks, blocks: blocks.length,
      uncertainBytes, uncertainRanges: ranges, sha256: sha },
  };
}

// --------------------------------------------------------------------------- channel simulation

/** Simulate synthesis + sequencing. Each designed strand is lost with probability ``loss``
 * (dropout), otherwise read ``reads`` times; every read gets independent per-base
 * substitutions/insertions/deletions. Output is shuffled, like real sequencing output. */
export function damage(text, { sub = 0, ins = 0, del = 0, loss = 0, reads = 1 }, rand = Math.random) {
  const stats = { strands: 0, lost: 0, reads: 0, substitutions: 0, insertions: 0, deletions: 0 };
  const out = [];
  for (const line of text.split('\n')) {
    if (!/^[ACGT]+$/.test(line)) continue;
    stats.strands++;
    if (rand() < loss) { stats.lost++; continue; }
    for (let c = 0; c < reads; c++) {
      // jump from error to error (geometric gaps) instead of drawing per base
      let s = '', pos = 0;
      const rate = sub + ins + del;
      while (rate > 0) {
        const gap = Math.floor(Math.log(1 - rand()) / Math.log(1 - Math.min(rate, 0.999999)));
        if (pos + gap >= line.length) break;
        s += line.slice(pos, pos + gap);
        pos += gap;
        const b = line[pos], r = rand() * rate;
        if (r < del) stats.deletions++;
        else if (r < del + ins) { stats.insertions++; s += BASES[Math.floor(rand() * 4)] + b; }
        else { stats.substitutions++; s += 'ACGT'.replace(b, '')[Math.floor(rand() * 3)]; }
        pos++;
      }
      s += line.slice(pos);
      out.push(s);
      stats.reads++;
    }
  }
  for (let i = out.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [out[i], out[j]] = [out[j], out[i]]; }
  return { text: out.map((s, i) => `>read_${i}\n${s}\n`).join(''), stats };
}
