#!/usr/bin/env python3
"""Bounded, additive record-scope diagnostic; never edits the frozen inputs."""
import argparse, csv, gzip, hashlib, json, os, subprocess, time
from pathlib import Path
import task12_scope_s1_s2 as old
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'results/ani-cases/2026-10-02/task12-astra-focused-scope'
BIN=ROOT/'.work/envs/task12-pyani/bin'
TARGETS={'GCA_001750725.1':('CP017280.1','GCA_019047005.1'), 'GCA_013731155.1':('CP056156.1','GCA_013743295.1')}
COMPARATORS={
'GCA_001750725.1':['GCA_019047005.1','GCA_000025565.1','GCA_009660455.1','GCA_001521715.1','GCA_001729765.1','GCA_001729805.1'],
'GCA_013731155.1':['GCA_013743295.1','GCA_001729705.1','GCA_001729725.1','GCA_001729785.1','GCA_001729745.1','GCA_040096145.1','GCA_052922205.1','GCA_019048245.1']}
rows=old.rows
write=old.write
sha=old.hashfile
fasta=old.fasta

def additive(original, whole, selected):
    assert selected in whole and selected not in original
    for rid,(_,seq) in original.items():
        assert rid in whole and seq==whole[rid][1],rid
    return {**original, selected:whole[selected]}

def write_fasta(p, records):
    with open(p,'w') as f:
        for rid,(h,seq) in records.items():
            f.write('>'+h+'\n')
            for start in range(0,len(seq),80):f.write(seq[start:start+80]+'\n')

def run_logged(argv, directory, env=None):
    directory.mkdir(parents=True,exist_ok=False)
    t=time.time()
    with open(directory/'stdout.log','w') as out,open(directory/'stderr.log','w') as err:
        p=subprocess.run(['/usr/bin/time','-v','-o',str(directory/'resources.txt')]+list(map(str,argv)),stdout=out,stderr=err,env=env)
    (directory/'command.json').write_text(json.dumps(dict(argv=list(map(str,argv)),exit_code=p.returncode,seconds=time.time()-t,job_id=os.getenv('SLURM_JOB_ID'),cpus=os.getenv('SLURM_CPUS_PER_TASK')),indent=2))
    assert p.returncode==0,(argv,p.returncode)

def prepare():
    tips={r['accession']:r for r in rows(old.AUD/'tip-scope-accounting.tsv')}
    manifest=[];batch=[];jobs=[]
    for a,(rid,neighbour) in TARGETS.items():
        c=ROOT/tips[a]['chrom_path'];w=ROOT/tips[a]['whole_path']
        original=fasta(c);whole=fasta(w);derived=additive(original,whole,rid)
        p=RUN/'inputs'/f'{a}_plus_{rid}.fna';assert not p.exists();write_fasta(p,derived)
        assert fasta(p)==derived
        d=dict(accession=a,flagged_record=rid,neighbour=neighbour,original_path=str(c),whole_path=str(w),derived_path=str(p),original_sha256=sha(c),whole_sha256=sha(w),derived_sha256=sha(p),original_bp=sum(len(v[1]) for v in original.values()),added_bp=len(whole[rid][1]),derived_bp=sum(len(v[1]) for v in derived.values()))
        manifest.append(d)
        # Frozen controls must be replayed to distinguish scope effects from predictor/training drift.
        batch.extend([(str(c),a+'_original'),(str(p),a+'_plus')])
        for b in COMPARATORS[a]:
            q=ROOT/tips[b]['chrom_path'];existing=None
            candidates=sorted(list((ROOT/'results/ani-cases').glob('*/*/anim/pairs/'+a+'__'+b+'/command.json'))+list((ROOT/'results/ani-cases').glob('*/*/anim/pairs/'+b+'__'+a+'/command.json')))
            for candidate in candidates:
                z=json.loads(candidate.read_text());args=z.get('argv',[])
                if z.get('exit_code')!=0 or '-i' not in args:continue
                src=Path(args[args.index('-i')+1]);fa=src/(a+'.fna');fb=src/(b+'.fna')
                if fa.exists() and fb.exists() and sha(fa)==sha(c) and sha(fb)==sha(q):existing=candidate.parent;break
            for mode,ap in [('original',c),('plus',p)]:
                reuse=str(existing) if mode=='original' and existing else ''
                jobs.append(dict(pair=a+'__'+b,mode=mode,target=a,comparator=b,target_path=str(ap),comparator_path=str(q),target_sha256=sha(ap),comparator_sha256=sha(q),reuse_source=reuse))
    write(RUN/'inputs.tsv',manifest)
    write(RUN/'anim-jobs.tsv',jobs)
    (RUN/'inputs/marker.batch.tsv').write_text(''.join(p+'\t'+a+'\n' for p,a in batch))
    # Archive baseline skani for predeclared comparators; also rank all neighbours for later inspection.
    census=[]
    for r in rows(old.PAIRS):
        if r['status']!='measured':continue
        for a,b,af,bf in [(r['a'],r['b'],r['af_a'],r['af_b']),(r['b'],r['a'],r['af_b'],r['af_a'])]:
            if a in TARGETS:census.append(dict(target=a,comparator=b,ani=r['ani'],target_af=af,comparator_af=bf,selected=b in COMPARATORS[a]))
    write(RUN/'original-skani-census.tsv',sorted(census,key=lambda x:(x['target'],-float(x['ani']),x['comparator'])))
    print(json.dumps(dict(pairs=len(jobs)//2,new_jobs=sum(not j['reuse_source'] for j in jobs),reused=sum(bool(j['reuse_source']) for j in jobs))))

def anim(index):
    j=rows(RUN/'anim-jobs.tsv')[index];assert not j['reuse_source']
    d=RUN/'anim'/j['mode']/j['pair'];inp=RUN/'inputs'/j['mode']/j['pair'];inp.mkdir(parents=True,exist_ok=False)
    for a,p,h in [(j['target'],j['target_path'],j['target_sha256']),(j['comparator'],j['comparator_path'],j['comparator_sha256'])]:
        assert sha(p)==h;(inp/(a+'.fna')).symlink_to(p)
    env=os.environ.copy();env.update(PATH=str(BIN)+':'+env['PATH'],OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLCONFIGDIR=str(RUN/'logs/mpl'))
    run_logged([BIN/'average_nucleotide_identity.py','-i',inp,'-o',d/'work','-m','ANIm','--workers','1','--nucmer_exe',BIN/'nucmer','--filter_exe',BIN/'delta-filter','--nocompress'],d,env)

def submit():
    specs=[]
    for i,j in enumerate(rows(RUN/'anim-jobs.tsv')):
        if not j['reuse_source']:specs.append((f'focus-ani-{i:02d}',f'python {ROOT}/scripts/task12_focused_scope.py anim --index {i}',1,'2G','15m'))
    specs.append(('focus-marker',f'bash {ROOT}/scripts/run_marker_cohort.sh --batchfile {RUN}/inputs/marker.batch.tsv --run-dir {RUN}/marker --threads 4',4,'12G','45m'))
    for name,command,cpus,mem,wall in specs:
        d=RUN/'scheduler'/name;d.mkdir(exist_ok=False)
        argv=['ssubmit','--json',name,command,'--mem',mem,'--time',wall,'--output',str(d/'job.%j.out'),'--error',str(d/'job.%j.err'),'--',f'--cpus-per-task={cpus}']
        for mode,op in [('--dry-run','plan'),('--test-only','test'),('','submit')]:
            cmd=argv.copy()
            if mode:cmd.insert(1,mode)
            r=subprocess.run(cmd,capture_output=True,text=True)
            (d/(op+'.json')).write_text(json.dumps(dict(argv=cmd,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr),indent=2))
            x=json.loads(r.stdout);assert r.returncode==0 and x['ok'] and x['schema_version']==1 and x['operation']==op,(name,r.stderr)
        print(name,x['submission']['job_id'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','anim','submit']);p.add_argument('--index',type=int);a=p.parse_args()
    {'prepare':prepare,'anim':lambda:anim(a.index),'submit':submit}[a.stage]()
