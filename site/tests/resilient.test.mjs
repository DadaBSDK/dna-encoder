// node --test site/tests/
import test from 'node:test';
import assert from 'node:assert/strict';
import { encodeResilient, decodeResilient, damage, planBlocks, ResilientError, STRAND_NT } from '../resilient.js';

function rng(seed) { // mulberry32
  return () => { seed |= 0; seed = (seed + 0x6d2b79f5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}
function bytes(n, seed) { const r = rng(seed); return Uint8Array.from({ length: n }, () => Math.floor(r() * 256)); }
const enc = (data, parity = 250, name = 'f.bin') => encodeResilient(new Blob([data]), name, parity);
const strands = (text) => text.split('\n').filter((l) => l && l[0] !== '>');

test('noiseless round trip across sizes and parity levels', async () => {
  for (const n of [0, 1, 31, 32, 33, 1000, 8191, 70000]) {
    for (const parity of [0, 100, 250, 1000]) {
      const data = bytes(n, n + parity);
      const { text, result } = await enc(data, parity);
      assert.ok(strands(text).every((s) => s.length === STRAND_NT));
      const out = decodeResilient(text);
      assert.equal(out.verified, true, `${n}/${parity}`);
      assert.deepEqual(out.bytes, data);
      assert.equal(out.report.uncertainBytes, 0);
      assert.equal(result.strands, strands(text).length);
    }
  }
});

test('no long homopolymers even for all-zero files', async () => {
  const { text } = await enc(new Uint8Array(5000));
  const longest = Math.max(...strands(text).map((s) => Math.max(...s.match(/(A+|C+|G+|T+)/g).map((r) => r.length))));
  assert.ok(longest <= 12, `longest run ${longest}`);
});

test('any p lost strands per block are rebuilt exactly', async () => {
  const data = bytes(40000, 7);
  const { text } = await enc(data, 250);
  const recs = strands(text);
  const blocks = planBlocks(Math.ceil(data.length / 32), 250);
  const header = recs.length - blocks.reduce((s, [k, p]) => s + k + p, 0);
  const r = rng(3);
  const keep = new Set(recs.keys());
  let q = header;
  for (const [k, p] of blocks) {
    const idx = Array.from({ length: k + p }, (_, i) => q + i).sort(() => r() - 0.5).slice(0, p);
    idx.forEach((i) => keep.delete(i));
    q += k + p;
  }
  const out = decodeResilient(recs.filter((_, i) => keep.has(i)).join('\n'));
  assert.equal(out.verified, true);
  assert.deepEqual(out.bytes, data);
});

test('one substitution, insertion or deletion in EVERY strand is repaired with zero parity', async () => {
  const data = bytes(6000, 11);
  const { text } = await enc(data, 0);
  const r = rng(5);
  const hit = strands(text).map((s) => {
    const i = Math.floor(r() * s.length), kind = Math.floor(r() * 3);
    if (kind === 0) return s.slice(0, i) + 'ACGT'.replace(s[i], '')[Math.floor(r() * 3)] + s.slice(i + 1);
    if (kind === 1) return s.slice(0, i) + 'ACGT'[Math.floor(r() * 4)] + s.slice(i);
    return s.slice(0, i) + s.slice(i + 1);
  });
  const out = decodeResilient(hit.join('\n'));
  assert.equal(out.verified, true);
  assert.equal(out.report.repaired + out.report.intact + out.report.skipped, hit.length);
  assert.equal(out.report.rejected, 0);
});

test('shuffled, lowercase, duplicated, junk and N-containing input still decodes', async () => {
  const data = bytes(3000, 21);
  const { text } = await enc(data, 250);
  const r = rng(9);
  const recs = strands(text).flatMap((s) => (r() < 0.3 ? [s, s] : [s]));
  recs.sort(() => r() - 0.5);
  recs.push('ACGTACGT', 'NNNN', 'hello');
  recs[0] = recs[0].slice(0, 50) + 'N' + recs[0].slice(51);
  const out = decodeResilient(recs.map((s) => (r() < 0.5 ? s.toLowerCase() : s)).join('\n'));
  assert.equal(out.verified, true);
});

test('compromise output: bytes not marked uncertain are always exact', async () => {
  for (let seed = 1; seed <= 12; seed++) {
    const r = rng(seed);
    const data = bytes(Math.floor(r() * 30000) + 1, seed);
    const { text } = await enc(data, [0, 100, 250][seed % 3]);
    const rates = { sub: r() * 0.03, ins: r() * 0.01, del: r() * 0.01, loss: r() * 0.4 };
    const hurt = damage(text, rates, r).text;
    let out;
    try { out = decodeResilient(hurt); } catch (e) { assert.ok(e instanceof ResilientError); continue; }
    assert.equal(out.bytes.length, data.length);
    const unsure = new Uint8Array(data.length);
    for (const [a, b] of out.report.uncertainRanges) unsure.fill(1, a, b);
    for (let i = 0; i < data.length; i++) if (!unsure[i]) assert.equal(out.bytes[i], data[i], `seed ${seed} byte ${i}`);
    const exact = data.every((b, i) => out.bytes[i] === b);
    assert.equal(out.verified, exact);
    if (out.report.uncertainBytes === 0) assert.equal(out.verified, true, `seed ${seed}`);
  }
});

test('realistic channel (1% errors/base, 3 reads, 5% dropout, 25% parity) recovers exactly', async () => {
  let exact = 0;
  for (let seed = 1; seed <= 10; seed++) {
    const data = bytes(20000, seed);
    const { text } = await enc(data, 250);
    const hurt = damage(text, { sub: 0.005, ins: 0.0025, del: 0.0025, loss: 0.05, reads: 3 }, rng(seed * 7)).text;
    exact += decodeResilient(hurt).verified;
  }
  assert.equal(exact, 10);
});

test('single read at 2% errors/base degrades to partial output instead of failing', async () => {
  const data = bytes(20000, 5);
  const { text } = await enc(data, 250);
  const out = decodeResilient(damage(text, { sub: 0.01, ins: 0.005, del: 0.005, loss: 0.05 }, rng(2)).text);
  assert.equal(out.verified, false);
  assert.ok(out.report.uncertainBytes > 0 && out.report.uncertainBytes < data.length);
  assert.equal(out.bytes.length, data.length);
});

test('decoding 1 MB with 1% substitutions and 3 reads finishes quickly', async () => {
  const data = bytes(1 << 20, 99);
  const { text } = await enc(data, 250);
  const hurt = damage(text, { sub: 0.01, reads: 3 }, rng(1)).text;
  const t0 = performance.now();
  const out = decodeResilient(hurt);
  const ms = performance.now() - t0;
  assert.equal(out.verified, true);
  assert.ok(ms < 20000, `${ms.toFixed(0)} ms`);
  console.log(`1 MB x3 reads: ${out.report.repaired} reads repaired, decode ${ms.toFixed(0)} ms`);
});

test('rejects text that is not an archive', () => {
  assert.throws(() => decodeResilient('>x\nACGTACGT\n'), ResilientError);
});
