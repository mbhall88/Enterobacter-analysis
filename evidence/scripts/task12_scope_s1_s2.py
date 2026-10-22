#!/usr/bin/env python3
"""Bounded diagnostic excluded-record/draft-anchor screen; never changes source scope."""
import argparse,csv,hashlib,json,subprocess,time
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('${PROJECT_ROOT}')
RUN=ROOT/'results/ani-cases/2026-09-30/task12-astra-scope-s1-s2'
AUD=ROOT/'results/project-review/2026-09-29/task18-claude-pre-leah-opus55/checks/scope-accounting'
PAIRS=ROOT/'results/ani-cases/2026-09-23/task12-astra-evidence/skani/pairs-complete.tsv'
GTDB=ROOT/'results/gtdb-context/2026-09-18/task05-opus-membership-subtree/tables/gtdb-membership.tsv'
BIN=ROOT/'.work/envs/task12-pyani/bin'
def rows(p):
    with open(p) as f: return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs,cols=None):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=cols or list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def fasta(p):
    d={}; name=None
    with open(p) as f:
        for l in f:
            if l.startswith('>'):
                name=l[1:].split()[0];assert name not in d;d[name]=[l[1:].strip(),[]]
            else:d[name][1].append(l.strip().upper())
    return {k:(h,''.join(s)) for k,(h,s) in d.items()}
def label(h):
    h=h.lower()
    return 'plasmid' if 'plasmid' in h else 'chromosome' if 'chromosome' in h or 'complete genome' in h else 'other'
def union(iv):
    out=[]
    for x,y in sorted((min(x,y),max(x,y)) for x,y in iv):
        if not out or x>out[-1][1]+1:out.append([x,y])
        else:out[-1][1]=max(y,out[-1][1])
    return out
def bp(iv):return sum(y-x+1 for x,y in union(iv))
def intersection(a,b):
    a=union(a);b=union(b);return bp([(max(x,u),min(y,v)) for x,y in a for u,v in b if max(x,u)<=min(y,v)])
def deletion_bounds(af_pct,total_bp,removed_bp):
    assert 0<=removed_bp<total_bp
    aligned=af_pct/100*total_bp
    return max(0,(aligned-removed_bp)/(total_bp-removed_bp)*100),min(100,aligned/(total_bp-removed_bp)*100)
def diagnostic(c,p):
    if c>=50 and p<20:return 'chromosome_like'
    if p>=50 and c<20:return 'plasmid_like'
    if c>=20 and p>=20:return 'shared_or_mixed'
    return 'unresolved_low_or_intermediate_coverage'
def hashfile(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def prepare():
    tips={r['accession']:r for r in rows(AUD/'tip-scope-accounting.tsv')}
    excluded=rows(AUD/'excluded-records.tsv');drafts=rows(AUD/'draft-type-size-vs-species.tsv')
    species={}
    for r in rows(GTDB):
        a=r['accession_version'];sp=r['gtdb_species']
        if sp:
            assert a not in species or species[a]==sp
            species[a]=sp
    wanted={r['accession'] for r in excluded if int(r['bp'])>=50000}|{r['accession'] for r in drafts}
    # One explicitly chromosome-labelled retained record; draft sequences cannot be neighbours.
    complete={a for a,t in tips.items() if t['chrom_records']=='1' and t['kept_record_labels']=='chromosome_label'}
    candidates=defaultdict(list)
    with open(PAIRS) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['status']!='measured':continue
            ani=float(r['ani']);af=min(float(r['af_a']),float(r['af_b']))
            for a,b in [(r['a'],r['b']),(r['b'],r['a'])]:
                if a in wanted and b in complete and a!=b and species.get(a) and species.get(a)==species.get(b):
                    candidates[a].append((ani,af,b))
    nb=[];jobs=[];inputs={str(AUD/n) for n in ['tip-scope-accounting.tsv','excluded-records.tsv','draft-type-size-vs-species.tsv']}|{str(PAIRS),str(GTDB)}
    for a in sorted(wanted):
        opts=sorted(candidates[a],key=lambda v:(-v[0],-v[1],v[2])); b=opts[0][2] if opts else ''
        nb.append(dict(accession=a,gtdb_species=species.get(a,''),eligible_neighbours=len(opts),neighbour=b,ani=opts[0][0] if opts else '',min_af=opts[0][1] if opts else '',status='selected' if b else 'no_eligible_complete_same_species_neighbour'))
        if not b:continue
        for kind in ['S1','S2']:
            wantedids={r['record'] for r in excluded if r['accession']==a and int(r['bp'])>=50000} if kind=='S1' else (set(fasta(tips[a]['chrom_path'])) if a in {r['accession'] for r in drafts} else set())
            if not wantedids:continue
            source=tips[a]['whole_path'] if kind=='S1' else tips[a]['chrom_path'];f=fasta(source)
            q=RUN/'inputs'/f'{kind}_{a}.fna'
            with open(q,'w') as o:
                for rec in sorted(wantedids):h,s=f[rec];o.write('>'+h+'\n'+s+'\n')
            target=tips[b]['whole_path'];inputs.update([source,target,tips[a]['chrom_path']])
            jobs.append(dict(kind=kind,accession=a,neighbour=b,query=str(q),target=target,self_chromosome=tips[a]['chrom_path']))
    write(RUN/'neighbours.tsv',nb);write(RUN/'alignment-jobs.tsv',jobs)
    write(RUN/'input-hashes.tsv',[dict(path=p,sha256=hashfile(p)) for p in sorted(inputs)])
    print(json.dumps(dict(wanted=len(wanted),selected=sum(bool(r['neighbour']) for r in nb),jobs=len(jobs),s1=sum(j['kind']=='S1' for j in jobs),s2=sum(j['kind']=='S2' for j in jobs))))
def command(args,out,err,log):
    t=time.time()
    with open(out,'w') as o,open(err,'w') as e:r=subprocess.run([str(x) for x in args],stdout=o,stderr=e)
    log.append(dict(argv=[str(x) for x in args],exit_status=r.returncode,elapsed_seconds=time.time()-t,stdout=str(out),stderr=str(err)))
    if r.returncode:raise RuntimeError(args)
def align_one(j):
    d=RUN/'alignments'/f"{j['kind']}_{j['accession']}";d.mkdir(exist_ok=True);log=[]
    try:
        for tag,target in [('neighbour',j['target'])]+([('self_chromosome',j['self_chromosome'])] if j['kind']=='S1' else []):
            p=d/tag
            # maxmatch and identity-only filtering keep shared/chromosome/plasmid HSPs visible.
            command([BIN/'nucmer','--maxmatch','-p',p,target,j['query']],str(p)+'.stdout',str(p)+'.stderr',log)
            command([BIN/'delta-filter','-i','95',str(p)+'.delta'],str(p)+'.i95.delta',str(p)+'.filter.stderr',log)
            command([BIN/'show-coords','-T','-r','-l','-H',str(p)+'.i95.delta'],str(p)+'.coords',str(p)+'.coords.stderr',log)
    finally:(d/'commands.json').write_text(json.dumps(log,indent=2))
    return j['accession']
def coords(p):
    data=defaultdict(lambda:defaultdict(list))
    with open(p) as f:
        for l in f:
            x=l.strip().split('\t')
            if not x:continue
            assert len(x)==11,(p,x)
            if float(x[6])<95:continue
            s,e=sorted((int(x[2]),int(x[3])));data[x[10]][x[9]].append((s,e))
    return data
def summarize(s2_only=False):
    jobs=rows(RUN/'alignment-jobs.tsv');nbs={r['accession']:r for r in rows(RUN/'neighbours.tsv')};results={};s2=[]
    for j in jobs:
        if s2_only and j['kind']!='S2':continue
        d=RUN/'alignments'/f"{j['kind']}_{j['accession']}";q=fasta(j['query']);targets=fasta(j['target']);iv=coords(d/'neighbour.coords')
        selfiv=coords(d/'self_chromosome.coords') if j['kind']=='S1' else {}
        for rec,(h,seq) in q.items():
            cats=defaultdict(list)
            for t,ints in iv[rec].items():cats[label(targets[t][0])]+=ints
            n=len(seq);c=bp(cats['chromosome']);p=bp(cats['plasmid']);o=bp(cats['other']);both=intersection(cats['chromosome'],cats['plasmid'])
            assert max(c,p,o)<=n
            r=dict(accession=j['accession'],record=rec,bp=n,neighbour=j['neighbour'],chromosome_aligned_bp=c,plasmid_aligned_bp=p,other_aligned_bp=o,shared_chromosome_plasmid_bp=both,chromosome_pct=100*c/n,plasmid_pct=100*p/n,classification=diagnostic(100*c/n,100*p/n),self_chromosome_aligned_bp=bp([z for ints in selfiv.get(rec,{}).values() for z in ints]) if j['kind']=='S1' else '',status='screened')
            if j['kind']=='S1':results[(j['accession'],rec)]=r
            else:s2.append(r)
    s1=[]
    for x in ([] if s2_only else rows(AUD/'excluded-records.tsv')):
        k=(x['accession'],x['record'])
        if k in results:r=results[k]
        else:
            r={k:'' for k in next(iter(results.values()))};r.update(accession=x['accession'],record=x['record'],bp=x['bp'],classification='not_screened',status='below_50000bp_bound' if int(x['bp'])<50000 else 'no_eligible_complete_same_species_neighbour')
        s1.append(r)
    for d in rows(AUD/'draft-type-size-vs-species.tsv'):
        a=d['accession']
        if not nbs[a]['neighbour']:
            tips={r['accession']:r for r in rows(AUD/'tip-scope-accounting.tsv')}
            for rec,(h,seq) in fasta(tips[a]['chrom_path']).items():
                r={k:'' for k in s2[0]};r.update(accession=a,record=rec,bp=len(seq),classification='unresolved_no_eligible_neighbour',status='no_eligible_complete_same_species_neighbour');s2.append(r)
    if not s2_only:write(RUN/'s1-records.tsv',s1)
    write(RUN/'s2-contigs.tsv',s2)
    anchors=[]
    for d in rows(AUD/'draft-type-size-vs-species.tsv'):
        rs=[r for r in s2 if r['accession']==d['accession']]
        assert len(rs)==int(d['contigs']) and sum(int(r['bp']) for r in rs)==int(d['draft_bp'])
        anchors.append(dict(**d,neighbour=nbs[d['accession']]['neighbour'],status=nbs[d['accession']]['status'],plasmid_like_contigs=sum(r['classification']=='plasmid_like' for r in rs),plasmid_like_bp=sum(int(r['bp']) for r in rs if r['classification']=='plasmid_like'),chromosome_like_contigs=sum(r['classification']=='chromosome_like' for r in rs),unresolved_or_shared_contigs=sum(r['classification'] not in ['chromosome_like','plasmid_like'] for r in rs)))
    write(RUN/'s2-anchors.tsv',anchors)
    assert s2_only or len(s1)==2839
    from collections import Counter
    print(json.dumps(dict(s1_status=Counter(r['status'] for r in s1),s1_class=Counter(r['classification'] for r in s1),s2_class=Counter(r['classification'] for r in s2)),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','align','summarize','summarize-s2']);p.add_argument('--workers',type=int,default=12);a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='summarize':summarize()
    elif a.stage=='summarize-s2':summarize(True)
    else:
        with ThreadPoolExecutor(max_workers=a.workers) as e:
            for result in e.map(align_one,rows(RUN/'alignment-jobs.tsv')):print(result,flush=True)
