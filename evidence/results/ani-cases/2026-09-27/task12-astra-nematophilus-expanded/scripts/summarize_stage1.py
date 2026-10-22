import csv,json,itertools,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rs):
 with open(p,'w') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
g={x['accession']:x for x in read(r/'inputs/panel-all.tsv')};new={x['accession'] for x in g.values() if x['role'] not in ['existing_focal','shared_neighbour_clade']};observations={}
for filename in ['skani-forward.tsv','skani-reverse.tsv']:
 for x in read(r/'outputs'/filename):
  a,b=Path(x['Ref_file']).stem,Path(x['Query_file']).stem
  if a==b:continue
  k=tuple(sorted([a,b]));v={'ani':float(x['ANI']),a:float(x['Align_fraction_ref']),b:float(x['Align_fraction_query'])};observations.setdefault(k,[]).append(v)
assert len(observations)==324
pairs=[]
for (a,b),vs in sorted(observations.items()):
 assert len(vs)>=2
 pairs.append(dict(a=a,b=b,ani=min(v['ani'] for v in vs),af_a=min(v[a] for v in vs),af_b=min(v[b] for v in vs),source='new_reciprocal_skani',observations=len(vs)))
for x in read('results/ani-cases/2026-09-23/task12-astra-evidence/skani/pairs-complete.tsv'):
 if x['a'] in g and x['b'] in g and x['a'] not in new and x['b'] not in new:
  pairs.append(dict(a=x['a'],b=x['b'],ani=float(x['ani']),af_a=float(x['af_a']),af_b=float(x['af_b']),source='frozen_skani',observations='reused'))
assert len(pairs)==820 and {(x['a'],x['b']) for x in pairs}==set(itertools.combinations(sorted(g),2))
for x in pairs:x['pass95_80']=x['ani']>=95 and min(x['af_a'],x['af_b'])>=80
write(r/'outputs/skani-pairs.tsv',pairs)
markers=read(r/'markers/identify/gtdbtk.bac120.markers_summary.tsv');write(r/'outputs/marker-counts.tsv',[{k:v for k,v in x.items() if not k.endswith('_names')} for x in markers])
ix={(x['a'],x['b']):x for x in pairs};summary=[]
for q in sorted(new):
 m=next(x for x in markers if x['name']==q)
 row=dict(accession=q,role=g[q]['role'],primary=g[q]['primary'],unique_markers=m['number_unique_genes'],multiple_markers=m['number_multiple_genes'])
 for label,t in [('bugandensis','GCA_019046905.1'),('nematophilus','GCA_026344075.1'),('ECC219','GCA_043966805.1')]:
  x=ix[tuple(sorted([q,t]))];row[label+'_ani']=x['ani'];row[label+'_min_af']=min(x['af_a'],x['af_b']);row[label+'_passes']=x['pass95_80']
 summary.append(row)
write(r/'outputs/additional-summary.tsv',summary)
flagged=[x['accession'] for x in summary if x['primary']=='True' and int(x['unique_markers'])<110]
assert flagged==['GCA_034107645.1'],flagged
qc=[]
for x in g.values():
 y=dict(x);y['local_primary']=x['primary']=='True' and x['accession'] not in flagged;y['local_qc_note']='102_unique_markers_hold_from_strict_primary_retain_in_sensitivity' if x['accession'] in flagged else 'as_planned';qc.append(y)
write(r/'inputs/panel-qc.tsv',qc)
assert sum(x['local_primary'] for x in qc)==38
(r/'outputs/stage1-validation.json').write_text(json.dumps({'status':'PASS','new_unordered_pairs':324,'reused_pairs':496,'all_pairs':820,'markers':len(markers),'planned_primary_marker_flags':flagged,'strict_local_primary':38,'118C5':'both retained for sensitivity; identity unresolved','dddh':'not_run'},indent=2)+'\n')
print(json.dumps(summary,indent=2))
