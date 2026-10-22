#!/usr/bin/env python3
"""Run frozen Group B versus hormaechei ANIm, reusing completed type comparisons."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def rows(p):
    with open(p) as f: return list(csv.DictReader(f,delimiter='\t'))


def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def write(p, rs):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rs)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--prepare',action='store_true')
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--pair',help='Run exactly this a__b pair in its own allocation')
    args=ap.parse_args(); root=Path.cwd(); run=args.run.resolve()
    src=root/'results/ani-cases/2026-09-23/task12-astra-evidence'
    cohort=root/'results/ani-cases/2026-09-24/task12-astra-group-b-hormaechei/outputs/genome-comparisons.tsv'
    envbin=root/'.work/envs/task12-pyani/bin'
    if args.prepare:
        assert not (run/'anim/selected-pairs.tsv').exists()
        genomes={r['accession']:r for r in rows(src/'inputs/resolved/genomes.tsv')}
        selected=[];inputs=[]
        for acc in sorted({r['accession'] for r in rows(cohort)}|{'GCA_019048245.1'}):
            r=genomes[acc];assert sha(r['source'])==r['sha256'],acc
            inputs.append(r)
        write(run/'inputs/genomes.tsv',inputs)
        for r in rows(cohort):
            a,b=sorted([r['accession'],'GCA_019048245.1']);key=a+'__'+b
            previous=src/'anim/pairs'/key
            selected.append(dict(a=a,b=b,selection_reasons='frozen_tree_group_B_to_hormaechei_type',
                                 reused=previous.exists(),source=str(previous) if previous.exists() else 'new'))
            if previous.exists():
                assert json.loads((previous/'command.json').read_text())['exit_code']==0
                (run/'anim/pairs'/key).symlink_to(previous)
            else:
                folder=run/'anim/fasta'/key;folder.mkdir()
                for acc in [a,b]:(folder/(acc+'.fna')).symlink_to(genomes[acc]['source'])
        assert len(selected)==438 and sum(r['reused'] for r in selected)==3
        write(run/'anim/selected-pairs.tsv',selected)
        (run/'inputs/selection-provenance.json').write_text(json.dumps({
            'cohort_sha256':sha(cohort),'resolved_manifest_sha256':sha(src/'inputs/resolved/genomes.tsv'),
            'genomes_hash_checked':len(inputs),'pairs':438,'reused':3,'new':435},indent=2)+'\n')
        return
    assert 1<=args.workers<=int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
    env=os.environ.copy();env['PATH']=str(envbin)+':'+env['PATH']
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
               MPLCONFIGDIR=str(run/'logs/matplotlib-cache'))
    def work(r):
        if r['reused']=='True':return
        key=r['a']+'__'+r['b'];dest=run/'anim/pairs'/key;dest.mkdir(exist_ok=False)
        argv=[str(envbin/'average_nucleotide_identity.py'),'-i',str(run/'anim/fasta'/key),
              '-o',str(dest/'work'),'-m','ANIm','--workers','1',
              '--nucmer_exe',str(envbin/'nucmer'),'--filter_exe',str(envbin/'delta-filter'),'--nocompress']
        record=dict(argv=argv,cwd=str(root),start_utc=datetime.now(timezone.utc).isoformat(),
                    slurm_job_id=os.getenv('SLURM_JOB_ID'),slurm_cpus=os.getenv('SLURM_CPUS_PER_TASK'))
        start=time.monotonic()
        with (dest/'stdout.log').open('w') as stdout,(dest/'stderr.log').open('w') as stderr:
            result=subprocess.run(['/usr/bin/time','-v','-o',str(dest/'resource.time')]+argv,
                                  env=env,stdout=stdout,stderr=stderr)
        record.update(exit_code=result.returncode,wall_seconds=time.monotonic()-start)
        (dest/'command.json').write_text(json.dumps(record,indent=2)+'\n')
        if result.returncode:raise RuntimeError('ANIm failed: '+key)
        print(key,round(record['wall_seconds'],2),flush=True)
    selected=rows(run/'anim/selected-pairs.tsv')
    if args.pair:
        selected=[r for r in selected if r['a']+'__'+r['b']==args.pair]
        assert len(selected)==1 and selected[0]['reused']=='False'
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(work,selected))
    print('Completed requested pair set.',flush=True)


if __name__=='__main__':main()
