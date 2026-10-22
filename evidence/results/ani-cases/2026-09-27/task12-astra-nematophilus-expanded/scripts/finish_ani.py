import csv,json,hashlib,statistics,itertools
from pathlib import Path
r=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rs):
 with open(p,'w') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
assert json.loads((r/'outputs/stage1-validation.json').read_text())['status']=='PASS'
pairs=read(r/'anim/selected-pairs.tsv');dirs=read(r/'anim/directional-audit/directional-alignments.tsv');assert len(pairs)==18 and len(dirs)==36
out=[]
for p in pairs:
 ds=[d for d in dirs if d['a']==p['a'] and d['b']==p['b']];assert len(ds)==2
 cmd=json.loads((r/'anim/pairs'/(p['a']+'__'+p['b'])/'command.json').read_text());assert cmd['exit_code']==0
 for d in ds:assert hashlib.sha256(Path(d['source']).read_bytes()).hexdigest()==d['sha256']
 ani=min(float(d['identity']) for d in ds);af=min(float(d[k]) for d in ds for k in ['reference_coverage','query_coverage'])
 out.append(dict(a=p['a'],b=p['b'],query=p['query'],reference=p['reference'],reason=p['selection_reasons'],anim=ani,min_four_af=af,passes95_80=ani>=95 and af>=80))
write(r/'outputs/anim-comparisons.tsv',out)
print(json.dumps(out,indent=2))
# Population summaries exclude poor assembly and both identity-unresolved/marker-flagged 118C5 records.
g={x['accession']:x for x in read(r/'inputs/panel-qc.tsv') if x['local_primary']=='True'}
group=lambda q:'nematophilus_candidate' if q=='GCA_026344075.1' else g[q]['group']
sets={k:{q for q in g if group(q)==k} for k in ['bugandensis','nematophilus_candidate','boundary','context']}
ps=read(r/'outputs/skani-pairs.tsv');sums=[]
for a,b in itertools.combinations_with_replacement(sets,2):
 chosen=[x for x in ps if x['a'] in g and x['b'] in g and sorted([group(x['a']),group(x['b'])])==sorted([a,b])]
 if chosen:sums.append(dict(group_a=a,group_b=b,pairs=len(chosen),min_ani=min(float(x['ani']) for x in chosen),median_ani=statistics.median(float(x['ani']) for x in chosen),max_ani=max(float(x['ani']) for x in chosen),passes95_80=sum(x['pass95_80']=='True' for x in chosen)))
write(r/'outputs/group-ani-summary.tsv',sums)
(r/'outputs/ani-phase-validation.json').write_text(json.dumps({'status':'PASS','anim_pairs':18,'directions':36,'skani_unique_pairs':820,'strict_local_panel':len(g),'groups':{k:len(v) for k,v in sets.items()},'dddh':'excluded'},indent=2)+'\n')
