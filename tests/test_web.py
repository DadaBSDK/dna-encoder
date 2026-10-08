import json
import threading
import time
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pytest

from dnastore.web import make_server, validate_options


@pytest.fixture
def server(tmp_path):
    srv = make_server(port=0, root=tmp_path)
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    yield srv, f'http://127.0.0.1:{srv.server_port}'
    srv.shutdown()
    srv.server_close()
    srv.workbench.executor.shutdown(wait=True)
    thread.join()


def submit(server, data, opts, name='hello.txt'):
    srv, base = server
    req = Request(base + '/api/run', data=data, headers={
        'X-Workbench-Token': srv.workbench.token, 'X-Options': quote(json.dumps(opts)), 'X-Filename': quote(name)})
    with urlopen(req) as response:
        ident = json.load(response)['id']
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline:
        with urlopen(base + '/api/jobs/' + ident) as response:
            job = json.load(response)
        if job['status'] in ('done', 'error'):
            return job
        time.sleep(.05)
    pytest.fail('Job did not finish')


def download(server, job, label):
    item = next(a for a in job['artifacts'] if a['label'] == label)
    with urlopen(server[1] + item['url']) as response:
        return response.read()


@pytest.mark.parametrize('mode', ['stream', 'oligo'])
def test_web_roundtrip(server, mode):
    data = b'Hello DNA!\n' * 8
    job = submit(server, data, {'mode': mode, 'action': 'roundtrip'})
    assert job['status'] == 'done', job
    assert job['result']['exact_match'] is True
    assert job['result']['input_name'] == 'hello.txt'
    assert download(server, job, 'Recovered file') == data
    report = json.loads(download(server, job, 'Experiment report'))
    assert report['verified'] is True


@pytest.mark.parametrize('mode,label,name', [('stream', 'DNA archive', 'archive.dna'), ('oligo', 'Oligo pool', 'pool.fasta')])
def test_encode_then_decode_uploaded_archive(server, mode, label, name):
    data = bytes(range(128))
    enc = submit(server, data, {'mode': mode, 'action': 'encode'}, 'test.bin')
    assert enc['status'] == 'done', enc
    archive = download(server, enc, label)
    dec = submit(server, archive, {'mode': mode, 'action': 'decode'}, name)
    assert dec['status'] == 'done', dec
    assert dec['result']['verified']
    assert download(server, dec, 'Recovered file') == data


def test_invalid_archive_and_traversal(server):
    job = submit(server, b'not an archive', {'mode': 'stream', 'action': 'decode'}, '../../bad.dna')
    assert job['status'] == 'error'
    assert not job['artifacts']
    assert (server[0].workbench.root / job['id'] / 'input' / 'bad.dna').exists()
    with pytest.raises(HTTPError) as exc:
        urlopen(server[1] + f'/api/files/{job["id"]}/input')
    assert exc.value.code == 404


def test_simulation_and_failed_recovery(server):
    job = submit(server, b'noise trial' * 20, {'mode': 'oligo', 'action': 'roundtrip', 'simulate': True,
                                            'channel': {'mean_coverage': 0}})
    assert job['status'] == 'done', job
    assert job['result']['verified'] is False
    assert job['result']['channel']['reads'] == 0
    assert all(a['label'] != 'Recovered file' for a in job['artifacts'])


def test_local_session_protection_and_static(server):
    _, base = server
    with urlopen(base) as response:
        assert b'Your files. Written in DNA.' in response.read()
    for request in [Request(base + '/api/run', data=b'abc'), Request(base + '/api/session', headers={'Host': 'evil.example'})]:
        with pytest.raises(HTTPError) as exc:
            urlopen(request)
        assert exc.value.code == 403
    with urlopen(base + '/api/session') as response:
        assert json.load(response)['token']


@pytest.mark.parametrize('raw', [{'mode':'invalid'}, {'mode':'oligo', 'oligo_len':99},
    {'mode':'oligo','channel':{'p_sub':float('nan')}}, {'mode':'oligo','channel':{'p_sub':.8,'p_del':.4}},
    {'mode':'oligo','simulate':'yes'}, {'mode':'oligo','codec':'missing'}])
def test_reject_invalid_settings(raw):
    with pytest.raises(ValueError):
        validate_options(raw)
