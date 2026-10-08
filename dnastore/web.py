"""Local browser workbench. Standard-library HTTP server, no frontend build step."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
from pathlib import Path
import secrets
import threading
import time
from urllib.parse import unquote, urlsplit
import uuid

STATIC = Path(__file__).parent / 'web_assets'
CODECS = {
    'goldman': ('goldman', {}),
    'whiten': ('naive2bit', {'whiten': True}),
    'steering-a': ('steering', {'alphabet': 'A', 'P': 8}),
    'steering-b': ('steering', {'alphabet': 'B', 'P': 6}),
    'steering-c': ('steering', {'alphabet': 'C', 'P': 8}),
    'fountain': ('fountain', {}),
}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as fh:
        while chunk := fh.read(65536):
            digest.update(chunk)
    return digest.hexdigest()


def validate_options(raw):
    if not isinstance(raw, dict):
        raise ValueError('settings must be an object')
    mode, action = raw.get('mode', 'stream'), raw.get('action', 'roundtrip')
    if mode not in ('stream', 'oligo') or action not in ('roundtrip', 'encode', 'decode'):
        raise ValueError('unknown mode or action')
    result = {'mode': mode, 'action': action}
    if mode == 'stream':
        return result
    from .channel import ChannelParams
    codec = raw.get('codec', 'goldman')
    if codec not in CODECS:
        raise ValueError('unknown codec')
    strength = raw.get('strength', 'erasure')
    if strength not in ('erasure', 'repair'):
        raise ValueError('invalid decoder strength')
    result.update(codec=codec, strength=strength)
    for name, default, low, high in [('oligo_len', 200, 100, 1000), ('parity', 150, 0, 1000), ('seed', 42, 0, 2**32 - 1)]:
        value = raw.get(name, default)
        if not isinstance(value, int) or isinstance(value, bool) or not low <= value <= high:
            raise ValueError(f'{name} must be an integer between {low} and {high}')
        result[name] = value
    simulate = raw.get('simulate', False)
    if not isinstance(simulate, bool):
        raise ValueError('simulate must be true or false')
    channel = raw.get('channel', {})
    if not isinstance(channel, dict) or set(channel) - {'p_sub', 'p_ins', 'p_del', 'dropout', 'mean_coverage'}:
        raise ValueError('invalid channel settings')
    for value in channel.values():
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise ValueError('channel settings must be numbers')
    params = ChannelParams.from_dict(channel)
    result.update(simulate=simulate, channel={key: getattr(params, key) for key in ('p_sub', 'p_ins', 'p_del', 'dropout', 'mean_coverage')})
    return result


class Workbench:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.token = secrets.token_urlsafe(32)
        self.jobs = {}
        self.lock = threading.Lock()
        self.executor = ThreadPoolExecutor(max_workers=1)

    def update(self, job, **values):
        with self.lock:
            job.update(values)

    def snapshot(self, job):
        with self.lock:
            return json.loads(json.dumps(job))

    def run(self, job, source, options):
        start = time.monotonic()
        root = source.parent.parent
        artifacts = []
        result = {'settings': options, 'input_bytes': source.stat().st_size, 'input_name': source.name}

        def artifact(path, label):
            artifacts.append({'name': path.name, 'label': label, 'bytes': path.stat().st_size,
                              'url': f'/api/files/{job["id"]}/{path.name}'})

        def phase(text):
            self.update(job, status='running', phase=text)

        try:
            mode, action = options['mode'], options['action']
            recovered = None
            if mode == 'stream':
                from .stream import encode_file, decode_file, MAGIC
                archive = source
                if action != 'decode':
                    phase('Encoding file into DNA')
                    archive = root / 'encoded.dna'
                    result['metrics'] = encode_file(source, archive)
                    artifact(archive, 'DNA archive')
                if action != 'encode':
                    phase('Recovering file and verifying SHA-256')
                    # Filename is for extension only, never interpreted as a path.
                    with archive.open('rb') as fh:
                        header = fh.readline(16385)
                    if not header.startswith(MAGIC):
                        raise ValueError('Upload a DNASTORE .dna archive for streaming decode')
                    name = json.loads(header[len(MAGIC):]).get('name', '')
                    suffix = Path(str(name)).suffix
                    suffix = suffix if len(suffix) <= 12 and suffix[1:].isalnum() else '.bin'
                    recovered = root / ('recovered' + suffix)
                    result['metrics'] = decode_file(archive, recovered)
                    result['verified'] = True
                with archive.open('rb') as fh:
                    fh.readline(16385)
                    preview = fh.read(4096).decode('ascii', errors='replace')
                result['sequence'] = ''.join(c for c in preview.split('!')[0] if c in 'ACGT')[:480]
            else:
                import numpy as np
                from .config import load_config
                from .primers import load_primers
                from .encoder import encode_bytes, write_fasta
                from .decoder import decode_reads, read_fasta
                from .metrics import pool_metrics
                from .screening import ScreeningConfig, screen_pool
                from .channel import ChannelParams, simulate
                name, params = CODECS[options['codec']]
                cfg = load_config(overrides={'workers': 1, 'oligo_len': options['oligo_len'],
                    'global_seed': options['seed'], 'codec': {'name': name, 'params': params},
                    'ecc': {'parity_permille': options['parity']}})
                primers = load_primers(cfg['primers'])
                if action != 'decode':
                    phase('Encoding error-corrected oligos')
                    enc = encode_bytes(source.read_bytes(), primers, cfg)
                    result['metrics'] = pool_metrics(enc)
                    phase('Screening oligo constraints')
                    screens = screen_pool([s for _, s in enc.oligos], primers, ScreeningConfig.from_dict(cfg['screening']), False, 1)
                    result['screening'] = {'violating_oligos': sum(bool(s['violations']) for s in screens),
                                           'oligos': len(screens)}
                    result['research_mode'] = 'Violations are reported; pool export is permitted for experiments.'
                    archive = root / 'pool.fasta'
                    write_fasta(enc.oligos, str(archive))
                    artifact(archive, 'Oligo pool')
                    reads = [s for _, s in enc.oligos]
                    result['sequence'] = ''.join(reads)[:480]
                    if action == 'roundtrip' and options['simulate']:
                        phase('Simulating sequencing errors')
                        records, truth = simulate(reads, ChannelParams.from_dict(options['channel']), np.random.default_rng(options['seed']))
                        reads = [r.seq for r in records]
                        result['channel'] = {'reads': len(reads), 'dropped_oligos': int(truth.dropped.sum()),
                            'substitutions': sum(r.n_sub for r in records), 'insertions': sum(r.n_ins for r in records),
                            'deletions': sum(r.n_del for r in records)}
                        readpath = root / 'reads.fasta'
                        write_fasta(list(enumerate(reads)), str(readpath))
                        artifact(readpath, 'Simulated reads')
                else:
                    phase('Reading FASTA / FASTQ')
                    reads = read_fasta(str(source))
                    result['sequence'] = ''.join(reads[:4])[:480]
                if action != 'encode':
                    phase('Decoding and checking file integrity')
                    decoded = decode_reads(reads, primers, cfg['oligo_len'], workers=1, strength=options['strength'], grouping='index')
                    result['decoder'] = decoded.report
                    result['verified'] = bool(decoded.ok)
                    if decoded.ok:
                        recovered = root / ('recovered' + decoded.extension)
                        recovered.write_bytes(decoded.data)
                    else:
                        result['recovery_failure'] = 'Decoder could not recover a verified file. See the decoder report.'
            if recovered is not None:
                artifact(recovered, 'Recovered file')
                result['recovered_sha256'] = sha(recovered)
                if action == 'roundtrip':
                    result['original_sha256'] = sha(source)
                    result['exact_match'] = result['original_sha256'] == result['recovered_sha256']
                    result['verified'] = result['verified'] and result['exact_match']
            result['elapsed_seconds'] = round(time.monotonic() - start, 3)
            report = root / 'report.json'
            report.write_text(json.dumps(result, indent=2))
            artifact(report, 'Experiment report')
            self.update(job, status='done', phase='Complete', result=result, artifacts=artifacts)
        except Exception as exc:
            self.update(job, status='error', phase='Experiment failed', error=str(exc), artifacts=artifacts)


def make_server(host='127.0.0.1', port=8765, root='results/workbench'):
    if host not in ('127.0.0.1', 'localhost'):
        raise ValueError('Workbench is local-only; use 127.0.0.1 or localhost')
    work = Workbench(root)

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def valid_host(self):
            return self.headers.get('Host') in (f'127.0.0.1:{self.server.server_port}', f'localhost:{self.server.server_port}')

        def send_data(self, data, code=200, kind='application/json'):
            if not isinstance(data, bytes):
                data = json.dumps(data).encode()
            self.send_response(code)
            self.send_header('Content-Type', kind)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Cache-Control', 'no-store')
            self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' blob:; media-src 'self' blob:; frame-ancestors 'none'; object-src 'none'")
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            if not self.valid_host():
                return self.send_data({'error': 'Invalid host'}, 403)
            path = urlsplit(self.path).path
            if path == '/api/session':
                with work.lock:
                    jobs = list(work.jobs.values())
                return self.send_data({'token': work.token, 'jobs': [work.snapshot(j) for j in jobs]})
            if path.startswith('/api/jobs/'):
                job = work.jobs.get(path.rsplit('/', 1)[-1])
                return self.send_data(work.snapshot(job) if job else {'error': 'Run not found'}, 200 if job else 404)
            if path.startswith('/api/files/'):
                parts = path.split('/')
                if len(parts) != 5 or parts[3] not in work.jobs:
                    return self.send_data({'error': 'File not found'}, 404)
                job = work.snapshot(work.jobs[parts[3]])
                item = next((a for a in job.get('artifacts', []) if a['url'] == path), None)
                if not item:
                    return self.send_data({'error': 'File not found'}, 404)
                file = work.root / parts[3] / item['name']
                self.send_response(200)
                self.send_header('Content-Type', 'application/octet-stream')
                self.send_header('Content-Disposition', f'attachment; filename="{file.name}"')
                self.send_header('Content-Length', str(file.stat().st_size))
                self.send_header('X-Content-Type-Options', 'nosniff')
                self.end_headers()
                with file.open('rb') as fh:
                    while chunk := fh.read(65536):
                        self.wfile.write(chunk)
                return
            names = {'/': 'index.html', '/app.js': 'app.js', '/style.css': 'style.css'}
            if path in names:
                file = STATIC / names[path]
                return self.send_data(file.read_bytes(), kind=mimetypes.guess_type(file.name)[0] or 'text/plain')
            self.send_data({'error': 'Not found'}, 404)

        def do_POST(self):
            if not self.valid_host() or not secrets.compare_digest(self.headers.get('X-Workbench-Token', ''), work.token):
                self.close_connection = True
                return self.send_data({'error': 'Invalid local session'}, 403)
            if self.path != '/api/run':
                return self.send_data({'error': 'Not found'}, 404)
            root = None
            try:
                options = validate_options(json.loads(unquote(self.headers.get('X-Options', '{}'))))
                if self.headers.get('Transfer-Encoding'):
                    raise ValueError('Content-Length is required')
                length = int(self.headers.get('Content-Length', '-1'))
                if length < 0:
                    raise ValueError('Content-Length is required')
                filename = Path(unquote(self.headers.get('X-Filename', 'upload.bin')).replace('\\', '/')).name
                if not filename or filename in ('.', '..'):
                    raise ValueError('Invalid filename')
                # Keep upload separate from output names, preserving suffix for gzip readers.
                ident = uuid.uuid4().hex
                root = work.root / ident
                root.mkdir()
                (root / "input").mkdir()
                source = root / "input" / filename
                remaining = length
                self.connection.settimeout(120)
                with source.open('xb') as fh:
                    while remaining:
                        chunk = self.rfile.read(min(65536, remaining))
                        if not chunk:
                            raise ValueError('Upload interrupted')
                        fh.write(chunk)
                        remaining -= len(chunk)
                job = {'id': ident, 'name': filename, 'status': 'queued', 'phase': 'Waiting for worker', 'created': time.time(), 'options': options}
                with work.lock:
                    work.jobs[ident] = job
                work.executor.submit(work.run, job, source, options)
                self.send_data({'id': ident}, 202)
            except (ValueError, OSError, TypeError) as exc:
                if root is not None:
                    import shutil
                    shutil.rmtree(root, ignore_errors=True)
                self.close_connection = True
                self.send_data({'error': str(exc)}, 400)

    server = ThreadingHTTPServer((host, port), Handler)
    server.workbench = work
    return server


def serve(host='127.0.0.1', port=8765, root='results/workbench'):
    server = make_server(host, port, root)
    print(f'DNA Workbench: http://127.0.0.1:{server.server_port}', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        server.workbench.executor.shutdown(wait=True)
