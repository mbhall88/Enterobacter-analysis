#!/usr/bin/env python3
"""Compare exact marker inputs and prepare paired fixed-column context trees."""
import gzip,json
from pathlib import Path
import task12_focused_scope as s
R=s.RUN;H=s.ROOT/'results/full-cohort/2026-09-22/task11-opus-full/marker'
def readfa(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt') as f:
  d={};k=None
  for l in f:
   if l.startswith('>'):k=l[1:].split()[0];assert k not in d;d[k]=''
   else:d[k]+=l.strip()
 return d

def main():
 checks=[];summaries=[];diffs=[]
 ms={mode:readfa(R/'marker'/mode/'align/gtdbtk.bac120.user_msa.fasta.gz') for mode in ['full','canonical']}
 hist={mode:readfa(H/'alignments'/f'{mode}.user_msa.fasta.gz') for mode in ms}
 cm=s.rows(H/'filtered/column-map.tsv');keep=[int(x['original_column'])-1 for x in cm if x['decision']=='retained']
 for a in s.TARGETS:
  for mode in ms:
   original,plus=ms[mode][a+'_original'],ms[mode][a+'_plus']
   assert original==hist[mode][a],(a,mode,'historical replay drift')
   use=range(len(original)) if mode!='full' else keep
   x=''.join(original[i] for i in use);y=''.join(plus[i] for i in use)
   summaries.append(dict(accession=a,alignment='primary_fixed_occupancy' if mode=='full' else 'canonical',sites=len(x),historical_control_identical=True,changed_columns=sum(u!=v for u,v in zip(x,y)),gaps_restored=sum(u=='-' and v!='-' for u,v in zip(x,y)),new_gaps=sum(u!='-' and v=='-' for u,v in zip(x,y)),aa_substitutions=sum(u!=v and u!='-' and v!='-' for u,v in zip(x,y))))
   if mode=='full':
    for i,(u,v) in enumerate(zip(original,plus)):
     if u!=v:diffs.append(dict(accession=a,original_column=i+1,marker=cm[i]['marker_id'],primary_retained=cm[i]['decision']=='retained',original=u,plus=v))
 s.write(R/'marker-input-sensitivity.tsv',summaries);s.write(R/'marker-column-changes.tsv',diffs)
 # Document marker locations in added record. Non-BAC120 top hits omitted explicitly.
 markerids={x['marker_id'] for x in cm};hits=[]
 for a,(rid,_) in s.TARGETS.items():
  p=R/'marker/identify/identify/intermediate_results/marker_genes'/f'{a}_plus'
  for source in ['pfam','tigrfam']:
   for r in s.rows(p/f'{a}_plus_{source}_tophit.tsv'):
    gene=r['Gene Id']
    for hit in r['Top hits (Family id,e-value,bitscore)'].split(';'):
     family,evalue,score=hit.split(',')
     if family in markerids and gene.startswith(rid+'_'):hits.append(dict(accession=a,record=rid,gene=gene,marker=family,evalue=evalue,bitscore=score))
 s.write(R/'added-record-bac120-hits.tsv',hits)
 # Bounded context: all approved type anchors + all18 prior oharae core members + S1 neighbours.
 panel={r['selected_accession_version'] for r in s.rows(s.ROOT/'manifests/types.tsv') if r['selected_accession_version']}
 sub=s.rows(s.ROOT/'results/ani-cases/2026-09-29/task12-astra-xiangfangensis-decision/approved-assignments.tsv')
 oharae={r['accession'] for r in sub if r['project_lineage']=='oharae'};assert len(oharae)==18
 panel|=oharae|set(s.TARGETS)|{v[1] for v in s.TARGETS.values()}
 assert panel<=set(hist['full'])
 s.write(R/'tree-panel.tsv',[dict(accession=a,oharae_member=a in oharae,scope_target=a in s.TARGETS) for a in sorted(panel)])
 (R/'trees').mkdir(exist_ok=False)
 for mask in ['primary','canonical']:
  source=hist['full' if mask=='primary' else 'canonical']
  for scope in ['original','plus']:
   path=R/'trees'/f'{mask}-{scope}.faa'
   with open(path,'w') as f:
    for a in sorted(panel):
     seq=ms['full' if mask=='primary' else 'canonical'][a+'_'+scope] if a in s.TARGETS else source[a]
     if mask=='primary':seq=''.join(seq[i] for i in keep)
     f.write('>'+a+'\n'+seq+'\n')
  x=readfa(R/'trees'/f'{mask}-original.faa');y=readfa(R/'trees'/f'{mask}-plus.faa')
  changed=[a for a in x if x[a]!=y[a]];assert changed==['GCA_013731155.1']
 (R/'marker-validation.json').write_text(json.dumps(dict(historical_controls_identical=4,tree_tips=len(panel),oharae_tips=len(oharae),only_changed_tree_row='GCA_013731155.1',summaries=summaries),indent=2))
 print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
