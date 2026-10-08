import { encodeBlob, decodeBlob, MAX_HEADER } from './stream.js';

const $ = (id) => document.getElementById(id);
const fmt = new Intl.NumberFormat('en');
const state = { file: null, archive: null, result: null, header: null, record: '', rest: null, urls: [] };

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

// --------------------------------------------------------------------------- encode

function chooseFile(file) {
  state.file = file;
  $('file-title').textContent = file.name;
  $('file-subtitle').textContent = `${bytes(file.size)} · ${file.type || 'unknown type'}`;
  $('encode-btn').disabled = false;
  showError($('encode-error'), file.size > 512 * 1024 * 1024
    ? 'Large file: the archive is about 4× this size and is held by your browser. It may be slow or run out of memory.' : '');
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
  prog.hidden = false;
  $('encode-badge').textContent = 'ENCODING';
  showError($('encode-error'), '');
  try {
    const { archive, result } = await encodeBlob(state.file, state.file.name, (f) => { prog.value = f; });
    state.archive = archive;
    state.result = result;
    await splitArchive(archive);
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

// Archive = header line + first record + rest. The lab edits only the first record.
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
  const { result, archive } = state;
  $('encode-empty').hidden = true;
  $('encode-result').hidden = false;
  $('encode-metrics').replaceChildren(
    metric('ORIGINAL', bytes(result.bytes), `${fmt.format(result.bytes)} bytes`),
    metric('NUCLEOTIDES', fmt.format(result.nucleotides), '4 per byte'),
    metric('ARCHIVE', bytes(archive.size), 'with framing'),
  );
  colorBases($('sequence'), state.record.slice(0, 480));
  $('sequence').hidden = !state.record;
  $('encode-sha').textContent = result.sha256;
  const a = $('download-archive');
  a.href = url(archive);
  a.download = `${result.name}.dna`;
  $('lab-empty').hidden = !!state.record;
  $('lab-empty').textContent = state.record ? '' : 'This file is empty, so there is no DNA record to edit. Try another file.';
  $('lab-body').hidden = !state.record;
  $('record').value = state.record;
  $('lab-result').hidden = true;
  updateDiff();
}

// --------------------------------------------------------------------------- mutation lab

function updateDiff() {
  const edited = $('record').value;
  const orig = state.record;
  let changed = 0;
  for (let i = 0; i < Math.max(edited.length, orig.length); i++) changed += edited[i] !== orig[i];
  const lenNote = edited.length !== orig.length ? `, length ${edited.length} vs ${orig.length}` : '';
  $('record-diff').textContent = changed ? `${changed} position${changed > 1 ? 's' : ''} changed${lenNote}` : 'unchanged';
  $('record-diff').className = changed ? 'changed' : '';
}

function mutate() {
  const box = $('record');
  const seq = box.value;
  if (!seq) return;
  const i = Math.floor(Math.random() * seq.length);
  const options = 'ACGT'.replace(seq[i], '');
  box.value = seq.slice(0, i) + options[Math.floor(Math.random() * options.length)] + seq.slice(i + 1);
  box.focus();
  box.setSelectionRange(i, i + 1);
  updateDiff();
}

async function labDecode() {
  const edited = $('record').value;
  const archive = new Blob([state.header, edited + '\n', state.rest]);
  const out = $('lab-result');
  out.hidden = false;
  out.replaceChildren(el('p', { class: 'hint' }, 'Decoding…'));
  out.replaceChildren(...(await decodeReport(archive, edited === state.record)));
}

async function decodeReport(archive, expectOk) {
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
      nodes.push(el('p', { class: 'hint' },
        'This archive has no redundancy, so the decoder can tell the data changed but not how to undo it. ',
        'It refuses the file instead of handing back silently damaged bytes.'));
    }
    return nodes;
  }
}

// --------------------------------------------------------------------------- decode panel

async function decodeUpload(file) {
  $('decode-title').textContent = `${file.name} · ${bytes(file.size)}`;
  const out = $('decode-result');
  out.hidden = false;
  out.replaceChildren(el('p', { class: 'hint' }, 'Decoding and verifying…'));
  out.replaceChildren(...(await decodeReport(file)));
}

// --------------------------------------------------------------------------- wiring

function dropzone(zone, input, onFile) {
  input.addEventListener('change', () => input.files[0] && onFile(input.files[0]));
  for (const ev of ['dragenter', 'dragover']) zone.addEventListener(ev, (e) => { e.preventDefault(); zone.classList.add('drag'); });
  for (const ev of ['dragleave', 'drop']) zone.addEventListener(ev, () => zone.classList.remove('drag'));
  zone.addEventListener('drop', (e) => { e.preventDefault(); if (e.dataTransfer.files[0]) onFile(e.dataTransfer.files[0]); });
}

dropzone($('dropzone'), $('file'), chooseFile);
dropzone($('decode-drop'), $('decode-file'), decodeUpload);
$('encode-btn').addEventListener('click', encode);
$('sample').addEventListener('click', () => {
  chooseFile(new File(['Hello from DNA! This sentence is stored as A, C, G and T.\n'], 'hello.txt', { type: 'text/plain' }));
  encode();
});
$('record').addEventListener('input', () => {
  const box = $('record');
  const pos = box.selectionStart;
  const clean = box.value.toUpperCase().replace(/\s+/g, '');
  if (clean !== box.value) { box.value = clean; box.setSelectionRange(pos, pos); }
  updateDiff();
});
$('mutate').addEventListener('click', mutate);
$('reset').addEventListener('click', () => { $('record').value = state.record; $('lab-result').hidden = true; updateDiff(); });
$('lab-decode').addEventListener('click', labDecode);
