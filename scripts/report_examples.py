"""Deterministic worked examples for the project report; no scientific ranking claims.

Run: .venv/bin/python scripts/report_examples.py --out results/report_examples_20261008
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
import time

import numpy as np
from dnastore.config import load_config
from dnastore.encoder import encode_bytes, write_fasta
from dnastore.decoder import decode_reads
from dnastore.primers import load_primers
from dnastore.channel import ChannelParams, simulate
from dnastore.metrics import pool_metrics
from dnastore.screening import ScreeningConfig, screen_pool
from dnastore.oligo import HEADER_INDEX_SPACE

ROOT=Path(__file__).resolve().parent.parent


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=ROOT/'results/report_examples_20261008')
    args=ap.parse_args()
    out=args.out
    if out.exists() and any(out.iterdir()):
        raise SystemExit('Use a new or empty output directory')
    out.mkdir(parents=True,exist_ok=True)
    text=('DNA storage project demonstration.\n'
          'We encode file bytes into DNA oligos, then recover the original bytes.\n'
          'Integrity check: exact byte equality and SHA-256.\n'
          'Unicode survives too: café, नमस्ते.\n')
    inputs={'text':text.encode('utf-8'),'binary':np.random.default_rng(20261008).bytes(1024)}
    rows=[]; details={'root_seed':20261008,'inputs':{},'pools':{},'examples':[]}
    pools={}
    for kind,data in inputs.items():
        cfg=load_config(ROOT/'configs/default.yaml', {'workers':1, 'global_seed':42,
            'codec':{'name':'goldman' if kind=='text' else 'naive2bit',
                     'params':{} if kind=='text' else {'whiten':True}}})
        primers=load_primers(cfg['primers'])
        enc=encode_bytes(data,primers,cfg)
        ext='txt' if kind=='text' else 'bin'
        (out/f'{kind}_input.{ext}').write_bytes(data)
        write_fasta(enc.oligos,str(out/f'{kind}_pool.fasta'))
        screens=screen_pool([s for _,s in enc.oligos],primers,ScreeningConfig.from_dict(cfg['screening']),False,1)
        details['inputs'][kind]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
        details['pools'][kind]={'config':cfg,'metrics':pool_metrics(enc),
            'violating_oligos':sum(bool(s['violations']) for s in screens),
            'first_data_oligo':next(s for idx,s in enc.oligos if idx>=HEADER_INDEX_SPACE)}
        pools[kind]=(enc,cfg,primers)
    cases=[('E1','Text, noiseless','text','clean',None),
           ('E2','Binary, noiseless','binary','clean',None),
           ('E3','Binary, two data oligos missing','binary','drop2',None),
           ('E4','Binary, 3% IDS at depth 12','binary','ids',ChannelParams(p_sub=.01,p_ins=.01,p_del=.01,
                mean_coverage=12,coverage='lognormal',lognormal_sigma=.27,rc_frac=.5,seed=101)),
           ('E5','Binary, 9% IDS at depth 2','binary','ids',ChannelParams(p_sub=.03,p_ins=.03,p_del=.03,
                mean_coverage=2,coverage='lognormal',lognormal_sigma=.27,rc_frac=.5,seed=202))]
    for eid,label,kind,channel,params in cases:
        enc,cfg,primers=pools[kind]
        original=inputs[kind]
        removed=[]; events={'substitutions':0,'insertions':0,'deletions':0}
        if channel=='ids':
            reads,truth=simulate([s for _,s in enc.oligos],params,np.random.default_rng(params.seed))
            seqs=[r.seq for r in reads]
            events={k:sum(getattr(r,attr) for r in reads) for k,attr in
                    [('substitutions','n_sub'),('insertions','n_ins'),('deletions','n_del')]}
        else:
            removed=[idx for idx,_ in enc.oligos if idx>=HEADER_INDEX_SPACE][:2] if channel=='drop2' else []
            seqs=[s for idx,s in enc.oligos if idx not in removed]
        write_fasta(list(enumerate(seqs)),str(out/f'{eid}_reads.fasta'))
        for strength in ('erasure','repair'):
            start=time.perf_counter()
            dec=decode_reads(seqs,primers,200,workers=1,strength=strength)
            elapsed=time.perf_counter()-start
            equal=dec.ok and dec.data==original
            recovered_hash=hashlib.sha256(dec.data).hexdigest() if dec.data is not None else None
            if equal:
                (out/f'{eid}_{strength}_recovered.{"txt" if kind=="text" else "bin"}').write_bytes(dec.data)
            row={'example':eid,'description':label,'strength':strength,'input_bytes':len(original),
                 'stored_bytes':enc.info.stored_len,'oligos':len(enc.oligos),'reads':len(seqs),
                 'removed_oligos':len(removed),'substitutions':events['substitutions'],
                 'insertions':events['insertions'],'deletions':events['deletions'],
                 'recovered':equal,'valid_data_oligos':dec.report.get('data_oligos_valid',0),
                 'repaired_oligos':dec.report.get('repaired',0),'error':dec.report.get('error',''),
                 'decode_s':elapsed}
            rows.append(row)
            details['examples'].append({**row,'removed_indices':removed,'channel_params':asdict(params) if params else None,
                'input_sha256':details['inputs'][kind]['sha256'],'output_sha256':recovered_hash,'decoder_report':dec.report})
            print(json.dumps(row),flush=True)
    with (out/'results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    details['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted((ROOT/'dnastore').rglob('*.py'))}
    details['example_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (out/'examples.json').write_text(json.dumps(details,indent=2)+'\n')
    print('Saved examples:',out)


if __name__=='__main__':
    main()
