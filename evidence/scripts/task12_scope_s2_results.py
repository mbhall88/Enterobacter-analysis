#!/usr/bin/env python3
"""Reparse every original/adjusted ANIm raw direction and preserve each denominator."""
import csv,json
from pathlib import Path
from pyani.anim import parse_delta
import task12_scope_s1_s2 as s
R=s.RUN
jobs=s.rows(R/'s2-anim-jobs.tsv');raw=[];summaries={}
for j in jobs:
 d=Path(j['reuse_source'] if j['reused']=='True' else j['output']);cmd=json.loads((d/'command.json').read_text());assert cmd['exit_code']==0
 work=d/'work';filters=sorted((work/'nucmer_output').glob('*/*.filter'));assert len(filters)==2
 a,b=j['type_accession'],j['query_accession'];lengths={a:sum(len(v[1]) for v in s.fasta(j['type_path']).values()),b:sum(len(v[1]) for v in s.fasta(j['query_path']).values())};out=[]
 with open(work/'ANIm_percentage_identity.tab') as f:
  rr=list(csv.reader(f,delimiter='\t'));matrix={r[0]:dict(zip(rr[0][1:],map(float,r[1:]))) for r in rr[1:]}
 for p in filters:
  ref,q=p.stem.split('_vs_');lr,lq,ani,errors=parse_delta(p);assert abs(matrix[ref][q]-ani)<1e-10
  r=dict(pair=j['pair'],mode=j['mode'],reference=ref,query=q,ani_pct=100*ani,reference_aligned_bp=lr,query_aligned_bp=lq,reference_denominator_bp=lengths[ref],query_denominator_bp=lengths[q],reference_af_pct=100*lr/lengths[ref],query_af_pct=100*lq/lengths[q],similarity_errors=errors,raw_filter=str(p),sha256=s.hashfile(p),reused=j['reused']);raw.append(r);out.append(r)
 minani=min(r['ani_pct'] for r in out);typeaf=min(r['reference_af_pct'] if r['reference']==a else r['query_af_pct'] for r in out);qaf=min(r['reference_af_pct'] if r['reference']==b else r['query_af_pct'] for r in out)
 summaries[(j['pair'],j['mode'])]=dict(ani=minani,type_af=typeaf,query_af=qaf,passes=minani>=95 and typeaf>=80 and qaf>=80)
s.write(R/'s2-anim-raw-directions.tsv',raw)
results=[]
for key in sorted({j['pair'] for j in jobs}):
 o=summaries[(key,'original')];a=summaries[(key,'adjusted')]
 results.append(dict(pair=key,original_ani_min_pct=o['ani'],adjusted_ani_min_pct=a['ani'],original_type_af_min_pct=o['type_af'],adjusted_type_af_min_pct=a['type_af'],original_query_af_min_pct=o['query_af'],adjusted_query_af_min_pct=a['query_af'],original_passes_95_80=o['passes'],adjusted_passes_95_80=a['passes'],ani_change_pp=a['ani']-o['ani'],type_af_change_pp=a['type_af']-o['type_af'],query_af_change_pp=a['query_af']-o['query_af'],pass_changed=a['passes']!=o['passes']))
s.write(R/'s2-anim-sensitivity.tsv',results)
(R/'s2-anim-validation.json').write_text(json.dumps(dict(pairs=len(results),directions=len(raw),identity_matrices_verified=True,original_and_derived_denominators_recomputed=True,pass_changes=sum(r['pass_changed'] for r in results)),indent=2))
print(json.dumps(results,indent=2))
