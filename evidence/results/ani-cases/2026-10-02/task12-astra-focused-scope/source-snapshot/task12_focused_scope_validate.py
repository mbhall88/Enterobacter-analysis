#!/usr/bin/env python3
"""Independent raw-delta/sequence checks without pyANI's reader or interval package."""
import csv,gzip,hashlib,json
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results/ani-cases/2026-10-02/task12-astra-focused-scope'
def rows(p):return list(csv.DictReader(open(p),delimiter='\t'))
def covered(xs):
    total=0;end=0
    for start,stop in sorted((min(a,b),max(a,b)) for a,b in xs):
        total+=max(0,stop-max(start-1,end));end=max(stop,end)
    return total

def delta(p):
    ref,qry={},{};bases=errors=0;r=q=None
    for line in Path(p).read_text().splitlines():
        x=line.split()
        if not x:continue
        if x[0].startswith('>'):r,q=x[0][1:],x[1]
        elif len(x)==7:
            a,b,c,d,e,_,_=map(int,x);ref.setdefault(r,[]).append((a,b));qry.setdefault(q,[]).append((c,d));bases+=abs(b-a)+abs(d-c)+2;errors+=e
    return sum(map(covered,ref.values())),sum(map(covered,qry.values())),100*(1-2*errors/bases),errors

def main():
    checks=0
    for r in rows(R/'anim-raw-directions.tsv'):
        a,b,ani,e=delta(r['raw_path']);assert a==int(r['reference_aligned_bp']) and b==int(r['query_aligned_bp']);assert abs(ani-float(r['ani_pct']))<1e-10;assert e==int(r['errors']);checks+=4
        assert 0<=float(r['reference_af'])<=100 and 0<=float(r['query_af'])<=100;checks+=1
    for r in rows(R/'inputs.tsv'):
        for scope in ['original','whole','derived']:
            assert hashlib.sha256(Path(r[scope+'_path']).read_bytes()).hexdigest()==r[scope+'_sha256'];checks+=1
        assert int(r['original_bp'])+int(r['added_bp'])==int(r['derived_bp']);checks+=1
    for mode in ['primary','canonical']:
        p=R/'trees'/f'{mode}-original.faa';q=R/'trees'/f'{mode}-plus.faa'
        def read(p):
            d={};k=None
            for l in p.read_text().splitlines():
                if l.startswith('>'):k=l[1:];d[k]=''
                else:d[k]+=l
            return d
        a,b=read(p),read(q);assert set(a)==set(b);assert [k for k in a if a[k]!=b[k]]==['GCA_013731155.1'];checks+=2
    out=dict(passed=True,independent_assertions=checks,raw_directions=56,input_hashes_unchanged=True,scope_addition_bp_conserved=True)
    (R/'independent-validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
if __name__=='__main__':main()
