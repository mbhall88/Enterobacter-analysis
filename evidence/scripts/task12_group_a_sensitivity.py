#!/usr/bin/env python3
"""Assess strict-core attrition with explicit diagnostic exclusions and pairwise deletion."""
import argparse,json,csv
from pathlib import Path
import numpy as np
from Bio import SeqIO
from task12_group_a_populations import BASE,TYPES,rows,write,fasta,strict_core
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output;dest=out/'coverage-sensitivity';dest.mkdir(exist_ok=False)
rr=rows(out/'inputs/membership.tsv');ids=[r['accession'] for r in rr];record={r['accession']:r for r in rr}
seq={r.id:np.frombuffer(str(r.seq).upper().encode(),dtype='S1') for r in SeqIO.parse(BASE/'core-genome/ledger-subset/alignment.fna','fasta') if r.id in set(ids)}
a=np.stack([seq[q] for q in ids]);del seq
valid=np.isin(a,np.frombuffer(b'ACGT',dtype='S1')); summaries=[]
for name,selected in [('without_high_missing',[q for q in ids if float(record[q]['missing_fraction'])<=.05 or q in TYPES]),('intersection',[q for q in ids if record[q]['in_all_three_minimal_type_clades']=='True'])]:
 ix=[ids.index(q) for q in selected];b=a[ix];_,keep,var=strict_core(b);snps=b[:,keep][:,var];fasta(dest/(name+'.fasta'),selected,snps)
 summaries.append(dict(dataset=name,n=len(selected),excluded=';'.join(q for q in ids if q not in selected),strict_core_bases=int(keep.sum()),variable_columns=int(var.sum())))
# Complete pairwise deletion: retain >=95% ACGT occupancy globally; unknowns never count as matches/mismatches.
keep=valid.mean(0)>=.95;a=a[:,keep];valid=valid[:,keep];matches=np.zeros((len(ids),len(ids)),dtype=np.float64);callable=matches.copy()
for start in range(0,a.shape[1],32768):
 b=a[:,start:start+32768];v=valid[:,start:start+32768].astype('float32');callable+=v@v.T
 for base in [b'A',b'C',b'G',b'T']:
  x=(b==base).astype('float32');matches+=x@x.T
assert np.all(callable>0) and np.all(matches<=callable)
distance=(callable-matches)/callable;np.savez_compressed(dest/'pairwise-deletion.npz',ids=np.array(ids),distance=distance,callable=callable,mismatches=callable-matches)
write(dest/'pairwise-deletion.tsv',[dict(a=ids[i],b=ids[j],callable_bases=int(callable[i,j]),mismatches=int(callable[i,j]-matches[i,j]),divergence=distance[i,j]) for i in range(len(ids)) for j in range(i+1,len(ids))])
write(dest/'datasets.tsv',summaries);(dest/'summary.json').write_text(json.dumps({'occupancy95_columns':a.shape[1],'minimum_pairwise_callable':int(callable.min()),'datasets':summaries},indent=2))
print((dest/'summary.json').read_text())
