import { encodeBlob, decodeBlob, MAGIC, MAX_HEADER } from './stream.js';
import { encodeResilient } from './resilient.js';

const $ = (id) => document.getElementById(id);
const fmt = new Intl.NumberFormat('en');
const state = { file: null, mode: 'resilient', archive: null, text: '', result: null, header: null, record: '', rest: null, urls: [] };

function bytes(n) {
  const units = ['B', 'KB', 'MB', 'GB'];
  let i = 0;
  while (n >= 1024 && i < units.length - 1) { n /= 1024; i++; }
  return `${i ? n.toFixed(1) : n} ${units[i]}`;
}

function url(blob) {
  const u = URL.createObjectURL(blob);
  state.urls.push(u);
  return u;
}

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') node.className = v; else node.setAttribute(k, v);
  }
  for (const c of children) node.append(c);
  return node;
}

function showError(box, message) {
  box.textContent = message;
  box.hidden = !message;
}

function colorBases(container, seq) {
  container.replaceChildren(...Array.from(seq, (b) => el('span', { class: b.toLowerCase() }, b)));
}

function metric(label, value, small) {
  return el('div', { class: 'metric' }, el('label', {}, label), el('strong', {}, value), el('small', {}, small || ''));
}

// --------------------------------------------------------------------------- helix, ticker, composition

const COMP = { A: 'T', C: 'G', G: 'C', T: 'A' };
const COLORS = { A: '#2f8a5f', C: '#c98a1b', G: '#4a63b5', T: '#c0504d' };
const SVG = 'http://www.w3.org/2000/svg';
const helix = { seq: '', phase: 0, nodes: null, frame: 0 };

function randomBases(n) {
  return Array.from({ length: n }, () => 'ACGT'[Math.floor(Math.random() * 4)]).join('');
}

function buildHelix(seq) {
  const svg = $('helix');
  const rungs = 18;
  helix.seq = seq.slice(0, rungs).padEnd(rungs, 'A');
  const s1 = document.createElementNS(SVG, 'path'), s2 = document.createElementNS(SVG, 'path');
  s1.setAttribute('class', 'strand'); s2.setAttribute('class', 'strand');
  s1.setAttribute('stroke', '#1d3b33'); s2.setAttribute('stroke', '#24654f');
  const rows = Array.from(helix.seq, (b) => {
    const g = document.createElementNS(SVG, 'g');
    const l1 = document.createElementNS(SVG, 'line'), l2 = document.createElementNS(SVG, 'line');
    l1.setAttribute('class', 'rung'); l2.setAttribute('class', 'rung');
    l1.setAttribute('stroke', COLORS[b]); l2.setAttribute('stroke', COLORS[COMP[b]]);
    const t1 = document.createElementNS(SVG, 'text'), t2 = document.createElementNS(SVG, 'text');
    t1.textContent = b; t2.textContent = COMP[b];
    t1.setAttribute('fill', COLORS[b]); t2.setAttribute('fill', COLORS[COMP[b]]);
    g.append(l1, l2, t1, t2);
    return { g, l1, l2, t1, t2 };
  });
  svg.replaceChildren(s1, ...rows.map((r) => r.g), s2);
  helix.nodes = { s1, s2, rows };
  drawHelix();
}

function drawHelix() {
  const { s1, s2, rows } = helix.nodes;
  const cx = 160, amp = 92, top = 24, step = 22, turn = 0.52;
  const pts = (sign) => {
    let d = '';
    for (let y = 0; y <= rows.length * step; y += 4) {
      const x = cx + sign * amp * Math.sin(helix.phase + (y / step) * turn);
      d += `${d ? 'L' : 'M'}${x.toFixed(1)},${(top + y - step / 2).toFixed(1)}`;
    }
    return d;
  };
  s1.setAttribute('d', pts(1));
  s2.setAttribute('d', pts(-1));
  rows.forEach((r, i) => {
    const a = helix.phase + (i + 0.5) * turn;
    const x1 = cx + amp * Math.sin(a), x2 = cx - amp * Math.sin(a), y = top + i * step;
    const depth = 0.35 + 0.65 * (Math.cos(a) * 0.5 + 0.5);
    r.l1.setAttribute('x1', x1); r.l1.setAttribute('x2', cx); r.l2.setAttribute('x1', cx); r.l2.setAttribute('x2', x2);
    for (const l of [r.l1, r.l2]) { l.setAttribute('y1', y); l.setAttribute('y2', y); }
    r.t1.setAttribute('x', x1 + Math.sign(x1 - cx || 1) * 11); r.t1.setAttribute('y', y);
    r.t2.setAttribute('x', x2 + Math.sign(x2 - cx || -1) * 11); r.t2.setAttribute('y', y);
    r.g.setAttribute('opacity', depth.toFixed(2));
  });
}

function animateHelix() {
  if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) return;
  let last = performance.now();
  const tick = (now) => {
    helix.phase += (now - last) * 0.0009;
    last = now;
    drawHelix();
    helix.frame = requestAnimationFrame(tick);
  };
  helix.frame = requestAnimationFrame(tick);
}

// Decorations must never stop encoding/decoding (e.g. if a stale cached page lacks an element).
function decorate(fn) {
  try { fn(); } catch (e) { console.warn('decoration skipped:', e); }
}

// Two identical halves scrolled by -50% loop seamlessly and fill the band from the first frame.
function fillTicker(seq) {
  let half = seq || randomBases(240);
  while (half.length < 240) half += seq;
  half = half.slice(0, 240) + '   ';
  colorBases($('ticker'), half + half);
}

function composition(seq) {
  const counts = { A: 0, C: 0, G: 0, T: 0 };
  for (const b of seq) counts[b] += 1;
  const n = seq.length || 1;
  $('composition').replaceChildren(...'ACGT'.split('').map((b) => {
    const pct = (100 * counts[b]) / n;
    const s = el('span', { title: `${b}: ${pct.toFixed(1)}%` });
    s.style.width = `${pct}%`;
    s.style.background = COLORS[b];
    return s;
  }));
  return counts;
}

// --------------------------------------------------------------------------- background decoding

let worker = null, queue = Promise.resolve();
function inWorker(message) {
  const run = async () => {
    try {
      worker ??= new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
      return await new Promise((resolve, reject) => {
        worker.onmessage = ({ data }) => resolve(data);
        worker.onerror = (e) => { e.preventDefault(); worker = null; reject(e); };
        worker.postMessage(message);
      });
    } catch {
      // module workers unavailable: run on the main thread instead
      const { decodeResilient, damage, ResilientError } = await import('./resilient.js');
      try {
        if (message.op === 'damage') return { ok: true, ...damage(message.text, message.channel) };
        const out = decodeResilient(message.text);
        return { ok: true, bytes: out.bytes, verified: out.verified, meta: out.meta, report: out.report };
      } catch (e) { return { ok: false, error: e instanceof ResilientError ? e.message : String(e) }; }
    }
  };
  queue = queue.then(run, run);
  return queue;
}

// --------------------------------------------------------------------------- encode

function currentMode() { return document.querySelector('input[name=mode]:checked').value; }

function syncMode() {
  const resilient = currentMode() === 'resilient';
  $('parity').hidden = $('parity-label').hidden = !resilient;
  $('encode-note').textContent = resilient
    ? 'Each 32-byte chunk becomes a 160-base strand with its index and CRC-32. Parity strands rebuild lost or damaged ones.'
    : 'Each byte becomes four bases (00=A · 01=C · 10=G · 11=T). Any changed base makes the whole file fail its SHA-256 check.';
}

function chooseFile(file) {
  state.file = file;
  $('file-title').textContent = file.name;
  $('file-subtitle').textContent = `${bytes(file.size)} · ${file.type || 'unknown type'}`;
  $('encode-btn').disabled = false;
  showError($('encode-error'), '');
  const preview = $('media-preview');
  preview.replaceChildren();
  preview.hidden = true;
  const kind = (file.type || '').split('/')[0];
  if (['image', 'video', 'audio'].includes(kind) && file.size < 64 * 1024 * 1024) {
    const tag = kind === 'image' ? 'img' : kind;
    const media = el(tag, { src: url(file), alt: `Preview of ${file.name}` });
    if (tag !== 'img') media.controls = true;
    preview.append(media);
    preview.hidden = false;
  }
}

async function encode() {
  if (!state.file) return;
  const btn = $('encode-btn'), prog = $('encode-progress');
  btn.disabled = true;
  $('encode-badge').textContent = 'ENCODING';
  showError($('encode-error'), '');
  state.mode = currentMode();
  try {
    if (state.mode === 'resilient') {
      const { archive, text, result } = await encodeResilient(state.file, state.file.name, Number($('parity').value));
      Object.assign(state, { archive, text, result });
      state.record = text.split('\n').find((l) => l && l[0] !== '>') || '';
    } else {
      prog.hidden = false;
      const { archive, result } = await encodeBlob(state.file, state.file.name, (f) => { prog.value = f; });
      Object.assign(state, { archive, result, text: '' });
      await splitArchive(archive);
    }
    renderEncoded();
    $('encode-badge').textContent = 'ENCODED';
  } catch (e) {
    showError($('encode-error'), e.message);
    $('encode-badge').textContent = 'FAILED';
  } finally {
    btn.disabled = false;
    prog.hidden = true;
  }
}

// Stream archive = header line + first record + rest. The lab edits only the first record.
async function splitArchive(archive) {
  const head = new Uint8Array(await archive.slice(0, MAX_HEADER + 1100).arrayBuffer());
  const h = head.indexOf(10);
  const r = head.indexOf(10, h + 1);
  const firstIsData = head[h + 1] !== 33; // '!' means an empty file: footer only
  state.header = archive.slice(0, h + 1);
  state.record = firstIsData ? new TextDecoder().decode(head.subarray(h + 1, r)) : '';
  state.rest = archive.slice(firstIsData ? r + 1 : h + 1);
}

function renderEncoded() {
  const { result, archive, mode } = state;
  $('encode-empty').hidden = true;
  $('encode-result').hidden = false;
  const perByte = result.bytes ? (result.nucleotides / result.bytes).toFixed(2) : '–';
  $('encode-metrics').replaceChildren(
    metric('ORIGINAL', bytes(result.bytes), `${fmt.format(result.bytes)} bytes`),
    mode === 'resilient'
      ? metric('STRANDS', fmt.format(result.strands), `${fmt.format(result.dataStrands)} data · ${fmt.format(result.parityStrands)} parity · ${result.headerStrands} header`)
      : metric('NUCLEOTIDES', fmt.format(result.nucleotides), '4 per byte'),
    mode === 'resilient'
      ? metric('NUCLEOTIDES', fmt.format(result.nucleotides), `${perByte} per byte incl. redundancy`)
      : metric('ARCHIVE', bytes(archive.size), 'with framing'),
  );
  const sample = mode === 'resilient'
    ? state.text.split('\n').filter((l) => l && l[0] !== '>').slice(0, 40).join('') : state.record;
  colorBases($('sequence'), sample.slice(0, 480));
  $('sequence').hidden = $('composition').hidden = !sample;
  if (sample) {
    const counts = composition(sample);
    const gc = (100 * (counts.G + counts.C)) / sample.length;
    decorate(() => {
      buildHelix(sample);
      fillTicker(sample);
      $('helix-caption').textContent = `fig. 1: the first ${helix.seq.length} bases of ${result.name} · GC ${gc.toFixed(0)}%`;
    });
  }
  $('encode-sha').textContent = result.sha256;
  const a = $('download-archive');
  a.href = url(archive);
  a.download = mode === 'resilient' ? `${result.name}.fasta` : `${result.name}.dna`;
  a.textContent = mode === 'resilient' ? '↓ Download strands (.fasta)' : '↓ Download .dna archive';
  setupLab();
}

// --------------------------------------------------------------------------- mutation lab

const VIEW_CHARS = 300_000;

function setupLab() {
  const resilient = state.mode === 'resilient';
  const hasDna = resilient ? !!state.text : !!state.record;
  $('lab-empty').hidden = hasDna;
  $('lab-empty').textContent = hasDna ? '' : 'This file is empty, so there is no DNA to damage. Try another file.';
  $('lab-body').hidden = !hasDna;
  $('sim').hidden = !resilient;
  $('manual').open = !resilient;
  $('lab-result').hidden = true;
  if (resilient) {
    let view = state.text;
    if (view.length > VIEW_CHARS) view = view.slice(0, view.lastIndexOf('\n>', VIEW_CHARS) + 1);
    state.view = view;
    state.rest = state.text.slice(view.length);
    $('record-title').textContent = state.rest ? 'FIRST STRANDS (FASTA)' : 'ALL STRANDS (FASTA)';
    $('record-help').textContent = 'Change, insert or delete bases anywhere, or delete whole strands. ' +
      (state.rest ? 'Only the first strands are shown; the rest are decoded unchanged.' : '');
  } else {
    state.view = state.record;
    $('record-title').textContent = 'FIRST DNA RECORD';
    $('record-help').textContent = 'Exact stream mode has no redundancy: any change makes the whole file fail. Switch to error-tolerant encoding to see recovery.';
  }
  $('record').value = state.view;
  updateDiff();
}

function updateDiff() {
  const edited = $('record').value, orig = state.view || '';
  let changed = 0;
  for (let i = 0; i < Math.max(edited.length, orig.length); i++) changed += edited[i] !== orig[i];
  const lenNote = edited.length !== orig.length ? `, length ${fmt.format(edited.length)} vs ${fmt.format(orig.length)}` : '';
  $('record-diff').textContent = changed ? `${fmt.format(changed)} position${changed > 1 ? 's' : ''} differ${lenNote}` : 'unchanged';
  $('record-diff').className = changed ? 'changed' : '';
}

function mutate() {
  const box = $('record'), text = box.value;
  const starts = [];
  let pos = 0;
  for (const line of text.split('\n')) {
    if (/^[ACGT]+$/.test(line)) starts.push([pos, line.length]);
    pos += line.length + 1;
  }
  if (!starts.length) return;
  const [s, len] = starts[Math.floor(Math.random() * starts.length)];
  const i = s + Math.floor(Math.random() * len);
  const options = 'ACGT'.replace(text[i], '');
  box.value = text.slice(0, i) + options[Math.floor(Math.random() * 3)] + text.slice(i + 1);
  box.focus();
  box.setSelectionRange(i, i + 1);
  updateDiff();
}

async function labDecode() {
  const out = $('lab-result');
  out.hidden = false;
  out.replaceChildren(busy('Decoding…'));
  if (state.mode === 'resilient') {
    const res = await inWorker({ op: 'decode', text: $('record').value + (state.rest || '') });
    out.replaceChildren(...renderRecovery(res));
  } else {
    const edited = $('record').value;
    const archive = new Blob([state.header, edited + '\n', state.rest]);
    out.replaceChildren(...(await streamReport(archive, edited === state.record)));
  }
  out.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

const PRESETS = {
  gentle: { sub: 0.2, ins: 0.05, del: 0.05, loss: 2, reads: 2 },
  realistic: { sub: 0.5, ins: 0.25, del: 0.25, loss: 5, reads: 3 },
  harsh: { sub: 1, ins: 0.5, del: 0.5, loss: 10, reads: 1 },
};
const SLIDERS = ['sub', 'ins', 'del', 'loss', 'reads'];

function syncSliders() {
  for (const k of SLIDERS) $(`${k}-out`).textContent = k === 'reads' ? `×${$(k).value}` : `${$(k).value}%`;
}

async function simulate() {
  const out = $('lab-result');
  out.hidden = false;
  const btn = $('simulate');
  btn.disabled = true;
  try {
    const channel = { sub: $('sub').value / 100, ins: $('ins').value / 100, del: $('del').value / 100,
      loss: $('loss').value / 100, reads: Number($('reads').value) };
    out.replaceChildren(busy('Damaging the DNA…'));
    const hurt = await inWorker({ op: 'damage', text: state.text, channel });
    out.replaceChildren(busy(`Decoding ${fmt.format(hurt.stats.reads)} damaged reads…`));
    const res = await inWorker({ op: 'decode', text: hurt.text });
    out.replaceChildren(...renderRecovery(res, hurt.stats));
    out.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } finally {
    btn.disabled = false;
  }
}

function busy(text) {
  return el('div', { class: 'busy' }, el('span', { class: 'spinner' }), text);
}

// --------------------------------------------------------------------------- results

function stat(label, value, cls = '') {
  return el('div', { class: `stat ${cls}` }, el('strong', {}, value), el('span', {}, label));
}

function isImage(bytesArr, name) {
  const b = bytesArr;
  return (b[0] === 0x89 && b[1] === 0x50) || (b[0] === 0xff && b[1] === 0xd8) || (b[0] === 0x47 && b[1] === 0x49)
    || (b[0] === 0x52 && b[8] === 0x57) || /\.(png|jpe?g|gif|webp|bmp)$/i.test(name);
}

function damageMap(total, ranges) {
  const bar = el('div', { class: 'damage-map', role: 'img', 'aria-label': `${ranges.length} uncertain regions` });
  for (const [a, b] of ranges) {
    const s = el('span');
    s.style.left = `${(100 * a) / total}%`;
    s.style.width = `max(2px, ${(100 * (b - a)) / total}%)`;
    bar.append(s);
  }
  return bar;
}

function renderRecovery(res, channel) {
  if (!res.ok) {
    return [el('div', { class: 'outcome failed' }, `✕ Nothing recoverable: ${res.error}.`),
      el('p', { class: 'hint' }, 'The metadata strands are stored 8 times, so this only happens with extreme damage. Try more reads per strand or less dropout.')];
  }
  const { report: r, meta, verified } = res;
  const total = meta.bytes, unsure = r.uncertainBytes;
  const exactPct = total ? (100 * (total - unsure)) / total : 100;
  const nodes = [];
  if (verified) {
    nodes.push(el('div', { class: 'outcome' }, `✓ Recovered exactly: ${fmt.format(total)} bytes, SHA-256 matches the original.`));
  } else if (unsure) {
    nodes.push(el('div', { class: 'outcome partial' },
      `◐ Partially recovered: ${exactPct.toFixed(1)}% of the file is guaranteed exact. ${fmt.format(unsure)} bytes in ${r.uncertainRanges.length} region${r.uncertainRanges.length > 1 ? 's' : ''} could not be rebuilt and hold best guesses.`));
  } else {
    nodes.push(el('div', { class: 'outcome failed' }, '✕ Integrity check failed although every strand verified. This should be extremely rare; please report it.'));
  }
  if (channel) {
    nodes.push(el('p', { class: 'channel' },
      `Channel: ${fmt.format(channel.strands)} strands, ${fmt.format(channel.lost)} dropped, ${fmt.format(channel.reads)} reads with `
      + `${fmt.format(channel.substitutions)} substitutions, ${fmt.format(channel.insertions)} insertions, ${fmt.format(channel.deletions)} deletions.`));
  }
  const lostStrands = r.guessed + r.missing;
  nodes.push(el('div', { class: 'stats' },
    stat('reads', fmt.format(r.candidates)),
    stat('intact', fmt.format(r.intact), 'ok'),
    stat('repaired', fmt.format(r.repaired), 'ok'),
    stat('merged by vote', fmt.format(r.consensus), 'ok'),
    stat('rebuilt by parity', fmt.format(r.rebuiltByParity), 'ok'),
    stat('unrecoverable', fmt.format(lostStrands), lostStrands ? 'bad' : '')));
  if (unsure) {
    nodes.push(el('div', { class: 'map-head' }, el('span', {}, 'Damage map'), el('span', {}, `${fmt.format(unsure)} uncertain bytes`)));
    nodes.push(damageMap(total, r.uncertainRanges));
  }
  const name = (meta.name || 'recovered.bin').split(/[\\/]/).pop();
  const file = new Blob([res.bytes]);
  if (isImage(res.bytes, name)) {
    const img = el('img', { src: url(file), alt: `Recovered ${name}` });
    const box = el('div', { class: 'preview-out' }, img);
    img.addEventListener('error', () => box.replaceChildren(el('p', { class: 'hint' },
      verified ? 'The browser cannot display this image format.'
        : 'These bytes do not form a displayable image yet: the damaged regions include parts the image format needs, such as its header. Raise the redundancy or the reads per strand, or download the partial file to inspect it.')));
    nodes.push(box);
  } else if (total && total < 2_000_000) {
    const text = new TextDecoder('utf-8').decode(res.bytes.subarray(0, 1200));
    if (!/[\u0000-\u0008\u000e-\u001f]/.test(text)) nodes.push(el('pre', { class: 'preview-text' }, text));
  }
  nodes.push(el('div', { class: 'downloads' },
    el('a', { href: url(file), download: verified ? name : `partial-${name}` }, verified ? `↓ Download ${name}` : `↓ Download partial ${name}`)));
  return nodes;
}

async function streamReport(archive, expectOk) {
  try {
    const { file, result } = await decodeBlob(archive);
    const name = result.name.split(/[\\/]/).pop() || 'recovered.bin';
    return [
      el('div', { class: 'outcome' }, `✓ Verified: ${fmt.format(result.bytes)} bytes recovered and the SHA-256 matches.`),
      el('div', { class: 'checks' }, el('p', {}, 'SHA-256', el('code', {}, result.sha256))),
      el('div', { class: 'downloads' }, el('a', { href: url(file), download: name }, `↓ Download ${name}`)),
    ];
  } catch (e) {
    const nodes = [el('div', { class: 'outcome failed' }, `✕ Rejected: ${e.message}`)];
    if (e.detail && e.detail.expected) {
      nodes.push(el('div', { class: 'checks' },
        el('p', {}, 'Fingerprint stored in the archive', el('code', {}, e.detail.expected)),
        el('p', {}, 'Fingerprint of the decoded bytes', el('code', { class: 'bad' }, e.detail.actual))));
    }
    if (expectOk === false) {
      nodes.push(el('p', { class: 'hint' }, 'Exact stream mode has no redundancy, so one changed base loses the file. Error-tolerant encoding survives this.'));
    }
    return nodes;
  }
}

// --------------------------------------------------------------------------- decode panel

async function decodeUpload(file) {
  $('decode-title').textContent = `${file.name} · ${bytes(file.size)}`;
  const out = $('decode-result');
  out.hidden = false;
  out.replaceChildren(busy('Decoding and verifying…'));
  const head = new TextDecoder().decode(await file.slice(0, MAGIC.length).arrayBuffer());
  if (head === MAGIC) out.replaceChildren(...(await streamReport(file)));
  else out.replaceChildren(...renderRecovery(await inWorker({ op: 'decode', text: await file.text() })));
}

// --------------------------------------------------------------------------- wiring

function dropzone(zone, input, onFile) {
  input.addEventListener('change', () => input.files[0] && onFile(input.files[0]));
  for (const ev of ['dragenter', 'dragover']) zone.addEventListener(ev, (e) => { e.preventDefault(); zone.classList.add('drag'); });
  for (const ev of ['dragleave', 'drop']) zone.addEventListener(ev, () => zone.classList.remove('drag'));
  zone.addEventListener('drop', (e) => { e.preventDefault(); if (e.dataTransfer.files[0]) onFile(e.dataTransfer.files[0]); });
}

decorate(() => {
  buildHelix(randomBases(18));
  fillTicker('');
  animateHelix();
});

dropzone($('dropzone'), $('file'), chooseFile);
dropzone($('decode-drop'), $('decode-file'), decodeUpload);
for (const r of document.querySelectorAll('input[name=mode]')) r.addEventListener('change', syncMode);
syncMode();
$('encode-btn').addEventListener('click', encode);
$('sample').addEventListener('click', () => {
  chooseFile(new File(['Hello from DNA! This sentence is stored as A, C, G and T, with enough redundancy to survive a few broken bases.\n'], 'hello.txt', { type: 'text/plain' }));
  encode();
});
$('record').addEventListener('input', () => {
  const box = $('record');
  if (state.mode === 'stream') {
    const pos = box.selectionStart;
    const clean = box.value.toUpperCase().replace(/\s+/g, '');
    if (clean !== box.value) { box.value = clean; box.setSelectionRange(pos, pos); }
  }
  updateDiff();
});
$('mutate').addEventListener('click', mutate);
$('reset').addEventListener('click', () => { $('record').value = state.view; $('lab-result').hidden = true; updateDiff(); });
$('lab-decode').addEventListener('click', labDecode);
$('simulate').addEventListener('click', simulate);
for (const k of SLIDERS) $(k).addEventListener('input', syncSliders);
for (const b of document.querySelectorAll('[data-preset]')) {
  b.addEventListener('click', () => { for (const [k, v] of Object.entries(PRESETS[b.dataset.preset])) $(k).value = v; syncSliders(); });
}
syncSliders();
