// Browser/Node port of dnastore/stream.py (DNASTORE-STREAM-1).
// Archives are byte-identical to the Python CLI's, and decoding applies the same checks:
// framing, record lengths, base alphabet, byte count and SHA-256, before any output exists.

export const MAGIC = '>DNASTORE-STREAM-1 ';
export const MAPPING = '00=A,01=C,10=G,11=T';
export const CHUNK_BYTES = 64 * 1024;
export const LINE_BYTES = 256;
export const MAX_HEADER = 16384;

export class StreamError extends Error {}

const BASES = [65, 67, 71, 84]; // A C G T
const BASE_VALUE = new Int8Array(256).fill(-1);
BASES.forEach((b, i) => { BASE_VALUE[b] = i; });
const enc = new TextEncoder();

// --------------------------------------------------------------------------- SHA-256 (incremental)

const K = new Uint32Array([
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]);

export class Sha256 {
  constructor() {
    this.h = new Uint32Array([0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]);
    this.buf = new Uint8Array(64);
    this.bufLen = 0;
    this.total = 0;
    this.w = new Uint32Array(64);
  }

  _block(p, o) {
    const w = this.w, h = this.h;
    for (let i = 0; i < 16; i++) w[i] = (p[o + 4 * i] << 24) | (p[o + 4 * i + 1] << 16) | (p[o + 4 * i + 2] << 8) | p[o + 4 * i + 3];
    for (let i = 16; i < 64; i++) {
      const a = w[i - 15], b = w[i - 2];
      const s0 = ((a >>> 7) | (a << 25)) ^ ((a >>> 18) | (a << 14)) ^ (a >>> 3);
      const s1 = ((b >>> 17) | (b << 15)) ^ ((b >>> 19) | (b << 13)) ^ (b >>> 10);
      w[i] = (w[i - 16] + s0 + w[i - 7] + s1) | 0;
    }
    let [a, b, c, d, e, f, g, hh] = h;
    for (let i = 0; i < 64; i++) {
      const S1 = ((e >>> 6) | (e << 26)) ^ ((e >>> 11) | (e << 21)) ^ ((e >>> 25) | (e << 7));
      const t1 = (hh + S1 + ((e & f) ^ (~e & g)) + K[i] + w[i]) | 0;
      const S0 = ((a >>> 2) | (a << 30)) ^ ((a >>> 13) | (a << 19)) ^ ((a >>> 22) | (a << 10));
      const t2 = (S0 + ((a & b) ^ (a & c) ^ (b & c))) | 0;
      hh = g; g = f; f = e; e = (d + t1) | 0; d = c; c = b; b = a; a = (t1 + t2) | 0;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
  }

  update(data) {
    let i = 0;
    this.total += data.length;
    if (this.bufLen) {
      while (i < data.length && this.bufLen < 64) this.buf[this.bufLen++] = data[i++];
      if (this.bufLen < 64) return this;
      this._block(this.buf, 0);
      this.bufLen = 0;
    }
    for (; i + 64 <= data.length; i += 64) this._block(data, i);
    while (i < data.length) this.buf[this.bufLen++] = data[i++];
    return this;
  }

  hex() {
    const bits = this.total * 8;
    const pad = new Uint8Array(((this.bufLen < 56 ? 56 : 120) - this.bufLen) + 8);
    pad[0] = 0x80;
    const hi = Math.floor(bits / 2 ** 32), lo = bits >>> 0;
    const n = pad.length;
    pad[n - 8] = hi >>> 24; pad[n - 7] = hi >>> 16; pad[n - 6] = hi >>> 8; pad[n - 5] = hi;
    pad[n - 4] = lo >>> 24; pad[n - 3] = lo >>> 16; pad[n - 2] = lo >>> 8; pad[n - 1] = lo;
    const total = this.total;
    this.update(pad);
    this.total = total;
    return Array.from(this.h, (x) => x.toString(16).padStart(8, '0')).join('');
  }
}

// --------------------------------------------------------------------------- JSON like Python's json.dumps

// Python json.dumps(ensure_ascii=True): escape '"', '\\' and everything outside 0x20..0x7e.
function pyString(s) {
  let out = '"';
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i);
    const ch = s[i];
    if (ch === '"') out += '\\"';
    else if (ch === '\\') out += '\\\\';
    else if (ch === '\n') out += '\\n';
    else if (ch === '\r') out += '\\r';
    else if (ch === '\t') out += '\\t';
    else if (ch === '\b') out += '\\b';
    else if (ch === '\f') out += '\\f';
    else if (c < 0x20 || c > 0x7e) out += '\\u' + c.toString(16).padStart(4, '0');
    else out += ch;
  }
  return out + '"';
}

function pyValue(v) {
  if (typeof v === 'string') return pyString(v);
  if (typeof v === 'number' && Number.isSafeInteger(v)) return String(v);
  throw new StreamError('unsupported metadata value');
}

// Python json.dumps(obj) with default separators (', ', ': '), keys in the given order.
function pyDumps(obj, keys) {
  return '{' + keys.map((k) => `${pyString(k)}: ${pyValue(obj[k])}`).join(', ') + '}';
}

// --------------------------------------------------------------------------- encode

const ENCODE = Array.from({ length: 256 }, (_, v) => [6, 4, 2, 0].map((s) => BASES[(v >> s) & 3]));

/** Encode a Blob/File into a DNASTORE-STREAM-1 archive Blob. ``onProgress(fraction)`` optional. */
export async function encodeBlob(blob, name, onProgress) {
  const metadata = { name: String(name), mapping: MAPPING };
  const header = enc.encode(MAGIC + pyDumps(metadata, ['name', 'mapping']) + '\n');
  if (header.length > MAX_HEADER) throw new StreamError('filename metadata too long');
  const parts = [header];
  const digest = new Sha256();
  let size = 0;
  for (let off = 0; off < blob.size; off += CHUNK_BYTES) {
    const chunk = new Uint8Array(await blob.slice(off, off + CHUNK_BYTES).arrayBuffer());
    digest.update(chunk);
    size += chunk.length;
    const lines = Math.ceil(chunk.length / LINE_BYTES);
    const out = new Uint8Array(chunk.length * 4 + lines);
    let o = 0;
    for (let s = 0; s < chunk.length; s += LINE_BYTES) {
      const end = Math.min(s + LINE_BYTES, chunk.length);
      for (let i = s; i < end; i++) {
        const e = ENCODE[chunk[i]];
        out[o++] = e[0]; out[o++] = e[1]; out[o++] = e[2]; out[o++] = e[3];
      }
      out[o++] = 10;
    }
    parts.push(out);
    if (onProgress) onProgress(Math.min(1, (off + chunk.length) / blob.size));
  }
  const result = { ...metadata, bytes: size, nucleotides: size * 4, sha256: digest.hex() };
  parts.push(enc.encode('!' + pyDumps(result, Object.keys(result).sort()) + '\n'));
  return { archive: new Blob(parts, { type: 'text/plain' }), result };
}

// --------------------------------------------------------------------------- decode

class LineReader {
  constructor(blob) { this.blob = blob; this.pos = 0; this.buf = new Uint8Array(0); this.start = 0; }

  async _fill() {
    if (this.pos >= this.blob.size) return false;
    const next = new Uint8Array(await this.blob.slice(this.pos, this.pos + CHUNK_BYTES).arrayBuffer());
    this.pos += next.length;
    const rest = this.buf.subarray(this.start);
    const merged = new Uint8Array(rest.length + next.length);
    merged.set(rest); merged.set(next, rest.length);
    this.buf = merged; this.start = 0;
    return true;
  }

  // Like Python's readline(limit): up to and including '\n', at most ``limit`` bytes.
  async readline(limit) {
    for (;;) {
      const avail = this.buf.subarray(this.start);
      const nl = avail.indexOf(10);
      if (nl >= 0 && nl < limit) { this.start += nl + 1; return avail.slice(0, nl + 1); }
      if (avail.length >= limit) { this.start += limit; return avail.slice(0, limit); }
      if (!(await this._fill())) { this.start += avail.length; return avail.slice(); }
    }
  }

  async atEnd() { return this.start >= this.buf.length && this.pos >= this.blob.size; }
}

function startsWith(bytes, text) {
  const t = enc.encode(text);
  return bytes.length >= t.length && t.every((b, i) => bytes[i] === b);
}

function parseJson(bytes) {
  try { return JSON.parse(new TextDecoder('utf-8', { fatal: true }).decode(bytes)); } catch { throw new StreamError('invalid archive metadata JSON'); }
}

function sameRecord(a, b) {
  const ka = Object.keys(a), kb = Object.keys(b);
  return ka.length === kb.length && ka.every((k) => Object.prototype.hasOwnProperty.call(b, k) && a[k] === b[k]);
}

/** Decode and verify an archive Blob. Returns ``{file: Blob, result}`` or throws StreamError.
 * On failure, ``err.detail`` may carry the first bad record (``{line, reason}``). */
export async function decodeBlob(blob, onProgress) {
  const r = new LineReader(blob);
  const header = await r.readline(MAX_HEADER + 1);
  if (header.length > MAX_HEADER || header[header.length - 1] !== 10 || !startsWith(header, MAGIC)) {
    throw new StreamError('not a supported DNASTORE streaming archive');
  }
  const metadata = parseJson(header.subarray(enc.encode(MAGIC).length));
  if (!metadata || typeof metadata !== 'object' || Array.isArray(metadata) || typeof metadata.name !== 'string' || metadata.mapping !== MAPPING) {
    throw new StreamError('invalid streaming metadata');
  }
  const digest = new Sha256();
  const parts = [];
  let size = 0, lineNo = 1;
  for (;;) {
    const line = await r.readline(MAX_HEADER + 1);
    lineNo += 1;
    if (!line.length) throw new StreamError('truncated archive: missing integrity footer');
    if (line.length > MAX_HEADER || line[line.length - 1] !== 10) throw new StreamError('invalid or oversized archive record');
    if (line[0] === 33) { // '!'
      const footer = parseJson(line.subarray(1));
      const expected = { ...metadata, bytes: size, nucleotides: size * 4, sha256: digest.hex() };
      if (!footer || typeof footer !== 'object' || Array.isArray(footer) || !sameRecord(footer, expected)) {
        const err = new StreamError('archive integrity check failed (metadata, length or SHA-256)');
        err.detail = { expected: footer && footer.sha256, actual: expected.sha256 };
        throw err;
      }
      if (!(await r.atEnd())) throw new StreamError('unexpected data after integrity footer');
      return { file: new Blob(parts), result: expected };
    }
    const n = line.length - 1;
    if (!n || n > LINE_BYTES * 4 || n % 4) throw new StreamError(`invalid DNA record length (line ${lineNo})`);
    const chunk = new Uint8Array(n / 4);
    for (let i = 0, j = 0; i < n; i += 4, j++) {
      const a = BASE_VALUE[line[i]], b = BASE_VALUE[line[i + 1]], c = BASE_VALUE[line[i + 2]], d = BASE_VALUE[line[i + 3]];
      if ((a | b | c | d) < 0) throw new StreamError(`DNA records must contain only uppercase A, C, G, T (line ${lineNo})`);
      chunk[j] = (a << 6) | (b << 4) | (c << 2) | d;
    }
    parts.push(chunk);
    digest.update(chunk);
    size += chunk.length;
    if (onProgress && blob.size) onProgress(Math.min(1, r.pos / blob.size));
  }
}
