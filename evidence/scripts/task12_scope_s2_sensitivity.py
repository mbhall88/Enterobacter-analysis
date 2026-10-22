#!/usr/bin/env python3
"""Select original near-AF-boundary comparisons and derivative ANIm pairs for S2."""
import argparse,csv,json,os,re,subprocess,time
from pathlib import Path
import task12_scope_s1_s2 as s
R=s.RUN
SRC=s.ROOT/'results/ani-cases/2026-09-23/task12-astra-evidence'
T04=s.ROOT/'results/ani-screen/2026-09-17/task04-opus-type-anchored/skani'
def acc(p):return re.search(r'GC[AF]_\d+\.\d+',str(p)).group()
def prepare():
    anchors={r['accession']:r for r in s.rows(R/'s2-anchors.tsv')};contigs=s.rows(R/'s2-contigs.tsv');tips={r['accession']:r for r in s.rows(s.AUD/'tip-scope-accounting.tsv')}
    derived={};meta=[]
    for a,r in anchors.items():
        remove={x['record'] for x in contigs if x['accession']==a and x['classification']=='plasmid_like'}
        if not remove:continue
        original=s.fasta(tips[a]['chrom_path']);p=R/'inputs'/f'{a}_diagnostic_plasmid_like_removed.fna'
        with open(p,'w') as f:
            for k,(h,seq) in original.items():
                if k not in remove:f.write('>'+h+'\n'+seq+'\n')
        removed=sum(len(original[k][1]) for k in remove);kept=sum(len(v[1]) for k,v in original.items() if k not in remove)
        assert removed==int(r['plasmid_like_bp']) and removed+kept==int(r['draft_bp']) and kept>0
        derived[a]=p;meta.append(dict(accession=a,original_path=tips[a]['chrom_path'],original_sha256=s.hashfile(tips[a]['chrom_path']),original_bp=removed+kept,removed_records=';'.join(sorted(remove)),removed_bp=removed,derived_bp=kept,derived_path=str(p),derived_sha256=s.hashfile(p)))
    s.write(R/'s2-derived-inputs.tsv',meta,['accession','original_path','original_sha256','original_bp','removed_records','removed_bp','derived_bp','derived_path','derived_sha256'])
    # Frozen full-cohort nonself ledger: all pairs per anchor; no preselection by a new alignment.
    census=[];select=[]
    with open(s.PAIRS) as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['status']!='measured':continue
            for a,b,af,otheraf in [(row['a'],row['b'],row['af_a'],row['af_b']),(row['b'],row['a'],row['af_b'],row['af_a'])]:
                if a not in anchors:continue
                near=75<=float(af)<=85
                why=('query_outside_frozen_full_manifest' if b not in tips else 'no_eligible_neighbour' if not anchors[a]['neighbour'] else 'no_plasmid_like_contigs' if a not in derived else 'type_af_outside_75_85' if not near else 'selected')
                r=dict(type_accession=a,query_accession=b,original_skani_ani=row['ani'],original_skani_type_af=af,original_skani_query_af=otheraf,near_boundary=near,disposition=why)
                census.append(r)
                if why=='selected':select.append(r)
    # For <1% proposed removal, retain all pairs as a fixed-alignment bound and
    # predefine extreme/boundary representatives; do not claim exhaustive realignment.
    reasons={}
    for a in derived:
        pool=[r for r in select if r['type_accession']==a]
        fraction=int(anchors[a]['plasmid_like_bp'])/int(anchors[a]['draft_bp'])
        if fraction>=0.01:
            for r in pool:reasons[(a,r['query_accession'])]=['all_AF_window_pairs_large_removal']
        else:
            rules=[('closest_AF_below_80',[r for r in pool if float(r['original_skani_type_af'])<80],lambda r:abs(float(r['original_skani_type_af'])-80)),
                   ('closest_AF_at_or_above_80',[r for r in pool if float(r['original_skani_type_af'])>=80],lambda r:abs(float(r['original_skani_type_af'])-80)),
                   ('closest_ANI_to_95',pool,lambda r:abs(float(r['original_skani_ani'])-95)),
                   ('highest_ANI',pool,lambda r:-float(r['original_skani_ani']))]
            for why,pop,key in rules:
                if pop:
                    rr=min(pop,key=lambda r:(key(r),r['query_accession']));reasons.setdefault((a,rr['query_accession']),[]).append(why)
    select=[]
    for r in census:
        a=r['type_accession'];L=int(anchors[a]['draft_bp']);removed=int(anchors[a]['plasmid_like_bp']);af=float(r['original_skani_type_af'])
        r['removed_bp']=removed
        r['fixed_alignment_type_af_lower_pct'],r['fixed_alignment_type_af_upper_pct']=s.deletion_bounds(af,L,removed)
        r['selection_reason']=';'.join(reasons.get((a,r['query_accession']),[]))
        if r['disposition']=='selected':
            if r['selection_reason']:select.append(r)
            else:r['disposition']='bounded_denominator_sensitivity_only_small_removal'
    s.write(R/'s2-pair-selection-census.tsv',census)
    # Each type derivative against a frozen query. For two affected types, preserve the other original input.
    jobs=[]
    for n,r in enumerate(select):
        a,b=r['type_accession'],r['query_accession'];key=f'{a}__{b}';q=Path(tips[b]['chrom_path']);orig=Path(tips[a]['chrom_path'])
        # Reuse only byte-identical original inputs, including completed later Task12 runs.
        existing=None
        for candidate in sorted(list((s.ROOT/'results/ani-cases').glob('*/*/anim/pairs/'+a+'__'+b+'/command.json'))+list((s.ROOT/'results/ani-cases').glob('*/*/anim/pairs/'+b+'__'+a+'/command.json'))):
            z=json.loads(candidate.read_text())
            if z.get('exit_code')!=0:continue
            argv=z.get('argv',[])
            if '-i' not in argv:continue
            inp=Path(argv[argv.index('-i')+1])
            fa,fb=inp/(a+'.fna'),inp/(b+'.fna')
            if fa.exists() and fb.exists() and s.hashfile(fa)==s.hashfile(orig) and s.hashfile(fb)==s.hashfile(q):
                existing=candidate.parent;break
        for mode,apath in [('original',orig),('adjusted',derived[a])]:
            dest=R/'sensitivity'/mode/key
            reuse=mode=='original' and existing is not None
            jobs.append(dict(pair=key,mode=mode,type_accession=a,query_accession=b,type_path=str(apath),query_path=str(q),output=str(dest),reused=str(reuse),reuse_source=str(existing) if reuse else ''))
    s.write(R/'s2-anim-jobs.tsv',jobs,['pair','mode','type_accession','query_accession','type_path','query_path','output','reused','reuse_source'])
    print(json.dumps(dict(affected_anchors=len(derived),selected_pairs=len(select),jobs=sum(j['reused']=='False' for j in jobs),reused=sum(j['reused']=='True' for j in jobs))))
def run(pair,mode):
    j=next(j for j in s.rows(R/'s2-anim-jobs.tsv') if j['pair']==pair and j['mode']==mode);assert j['reused']=='False'
    d=Path(j['output']);d.mkdir(parents=True,exist_ok=False);inp=d/'fasta';inp.mkdir()
    for a,p in [(j['type_accession'],j['type_path']),(j['query_accession'],j['query_path'])]:(inp/(a+'.fna')).symlink_to(Path(p).resolve())
    env=os.environ.copy();env.update(PATH=str(s.BIN)+':'+env['PATH'],OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    cmd=[str(s.BIN/'average_nucleotide_identity.py'),'-i',str(inp),'-o',str(d/'work'),'-m','ANIm','--workers','1','--nucmer_exe',str(s.BIN/'nucmer'),'--filter_exe',str(s.BIN/'delta-filter'),'--nocompress'];t=time.time()
    with open(d/'stdout.log','w') as out,open(d/'stderr.log','w') as err:r=subprocess.run(['/usr/bin/time','-v','-o',str(d/'resources.txt')]+cmd,env=env,stdout=out,stderr=err)
    (d/'command.json').write_text(json.dumps(dict(argv=cmd,exit_code=r.returncode,elapsed_seconds=time.time()-t,job_id=os.getenv('SLURM_JOB_ID'),cpus=os.getenv('SLURM_CPUS_PER_TASK'),type_sha256=s.hashfile(j['type_path']),query_sha256=s.hashfile(j['query_path'])),indent=2));assert r.returncode==0
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run']);p.add_argument('--pair');p.add_argument('--mode');a=p.parse_args();prepare() if a.stage=='prepare' else run(a.pair,a.mode)
