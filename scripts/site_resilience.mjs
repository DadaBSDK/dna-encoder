// Benchmark of the browser error-tolerant codec (site/resilient.js, DNASTORE-RS-1).
//
//   node scripts/site_resilience.mjs [out_dir]
//
// For every channel cell (per-base error rate x reads per strand x parity) it encodes
// independent random payloads, applies site/resilient.js damage() (5% strand dropout,
// error split 50% substitution / 25% insertion / 25% deletion), decodes, and records exact
// recovery, guaranteed-exact bytes, actually-correct bytes and decode time. All randomness
// is seeded, so a rerun reproduces the same numbers (timings aside).
import { mkdirSync, writeFileSync, readFileSync, existsSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { cpus } from 'node:os';
import { encodeResilient, decodeResilient, damage, ResilientError } from '../site/resilient.js';

const OUT = process.argv[2] || 'results/site_resilience';
if (existsSync(`${OUT}/results.csv`)) throw new Error(`${OUT} already has results; choose a new directory`);
mkdirSync(OUT, { recursive: true });

const PAYLOAD_BYTES = 20000;
const SEEDS = 10;
const ERROR_RATES = [0.003, 0.01, 0.02];
const READS = [1, 2, 3];
const PARITY = [100, 250, 500];
const DROPOUT = 0.05;

function rng(seed) { // mulberry32
  return () => { seed |= 0; seed = (seed + 0x6d2b79f5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}

const rows = [];
const t0 = Date.now();
for (const parity of PARITY) {
  for (const err of ERROR_RATES) {
    for (const reads of READS) {
      for (let s = 0; s < SEEDS; s++) {
        const seed = 1000003 * (s + 1) + 7919 * reads + Math.round(err * 1e5) + parity;
        const r = rng(seed);
        const data = Uint8Array.from({ length: PAYLOAD_BYTES }, () => Math.floor(r() * 256));
        const { text, result } = await encodeResilient(new Blob([data]), 'payload.bin', parity);
        const channel = { sub: err * 0.5, ins: err * 0.25, del: err * 0.25, loss: DROPOUT, reads };
        const hurt = damage(text, channel, rng(seed ^ 0x5bd1e995));
        const t = performance.now();
        let out = null, error = '';
        try { out = decodeResilient(hurt.text); } catch (e) { if (!(e instanceof ResilientError)) throw e; error = e.message; }
        const ms = performance.now() - t;
        const correct = out ? data.reduce((c, b, i) => c + (out.bytes[i] === b), 0) : 0;
        if (out) {  // safety property: bytes not reported uncertain must be exact
          const unsure = new Uint8Array(data.length);
          for (const [a, b] of out.report.uncertainRanges) unsure.fill(1, a, b);
          for (let i = 0; i < data.length; i++) if (!unsure[i] && out.bytes[i] !== data[i]) throw new Error(`uncertainty violated at seed ${seed}`);
        }
        rows.push({
          parity_permille: parity, error_per_base: err, reads_per_strand: reads, dropout: DROPOUT, seed,
          payload_bytes: PAYLOAD_BYTES, strands: result.strands, nt_per_byte: (result.nucleotides / PAYLOAD_BYTES).toFixed(4),
          exact: out ? Number(out.verified) : 0, header_lost: Number(!out),
          guaranteed_exact_frac: out ? (1 - out.report.uncertainBytes / PAYLOAD_BYTES).toFixed(5) : '0',
          actually_correct_frac: (correct / PAYLOAD_BYTES).toFixed(5),
          reads_total: hurt.stats.reads, intact: out?.report.intact ?? 0, repaired: out?.report.repaired ?? 0,
          consensus: out?.report.consensus ?? 0, rebuilt_by_parity: out?.report.rebuiltByParity ?? 0,
          unrecoverable_strands: out ? out.report.guessed + out.report.missing : '', decode_ms: ms.toFixed(1), error,
        });
      }
      process.stdout.write(`parity ${parity} err ${err} reads ${reads}: ${rows.slice(-SEEDS).reduce((c, x) => c + x.exact, 0)}/${SEEDS} exact\n`);
    }
  }
}

const cols = Object.keys(rows[0]);
writeFileSync(`${OUT}/results.csv`, [cols.join(','), ...rows.map((x) => cols.map((c) => JSON.stringify(x[c] ?? '')).join(','))].join('\n') + '\n');
const summary = [];
for (const parity of PARITY) for (const err of ERROR_RATES) for (const reads of READS) {
  const g = rows.filter((x) => x.parity_permille === parity && x.error_per_base === err && x.reads_per_strand === reads);
  const mean = (k) => g.reduce((c, x) => c + Number(x[k]), 0) / g.length;
  summary.push({ parity_permille: parity, error_per_base: err, reads_per_strand: reads, trials: g.length,
    exact: g.reduce((c, x) => c + x.exact, 0), header_lost: g.reduce((c, x) => c + x.header_lost, 0),
    mean_guaranteed_exact: +mean('guaranteed_exact_frac').toFixed(5), mean_actually_correct: +mean('actually_correct_frac').toFixed(5),
    nt_per_byte: +g[0].nt_per_byte, median_decode_ms: +[...g.map((x) => +x.decode_ms)].sort((a, b) => a - b)[g.length >> 1].toFixed(1) });
}
const sc = Object.keys(summary[0]);
writeFileSync(`${OUT}/summary.csv`, [sc.join(','), ...summary.map((x) => sc.map((c) => x[c]).join(','))].join('\n') + '\n');
const sha = (p) => createHash('sha256').update(readFileSync(p)).digest('hex');
writeFileSync(`${OUT}/provenance.json`, JSON.stringify({
  script: 'scripts/site_resilience.mjs', codec: 'site/resilient.js', codec_sha256: sha('site/resilient.js'),
  script_sha256: sha('scripts/site_resilience.mjs'), node: process.version, cpu: cpus()[0]?.model, started_unix_ms: t0,
  wall_s: (Date.now() - t0) / 1000, design: { PAYLOAD_BYTES, SEEDS, ERROR_RATES, READS, PARITY, DROPOUT,
    error_split: 'sub 50% / ins 25% / del 25% of the per-base rate; independent per read; strand dropout before reading' },
}, null, 2) + '\n');
console.log(`done: ${rows.length} trials in ${((Date.now() - t0) / 1000).toFixed(0)} s -> ${OUT}`);
