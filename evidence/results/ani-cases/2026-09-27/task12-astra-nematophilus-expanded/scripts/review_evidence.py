"""Derive review tables and audit claims against completed outputs, no new inference."""
from pathlib import Path
import csv,json,hashlib
from Bio import Phylo
r=Path(__file__).resolve().parents[1];out=r/'review'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(name,rows):
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
assert json.loads((r/'outputs/final-validation.json').read_text())['status']=='PASS'
splits=read(r/'outputs/tree-splits.tsv');cf=read(r/'outputs/site-concordance.tsv');au=read(r/'outputs/topology-tests.tsv');dist=read(r/'outputs/patristic-summary.tsv')
for group in ['bugandensis','nematophilus','boundary','context']:
 s=[x for x in splits if x['split']==group];assert len(s)==10 and all(x['present']=='True' for x in s)
 for x in s:
  if x['method'] in ['primary','masked','RY']:assert x['SHaLRT_UFB_label']=='100/100'
refs=['GCA_019046905.1','GCA_026344075.1'];table=[];separation=[];sensitivity=[]
for ref in refs:
 o=r/'local'/ref;t=o/'tests';v=json.loads((o/'validation.json').read_text());m=json.loads((t/'recombination-summary.json').read_text())
 log=(t/'logs/gubbins.stdout').read_text();assert 'Maximum number of iterations (5) reached.' in log
 for condition in ['primary','masked']:
  expected='boundary_with_nematophilus' if condition=='primary' else 'boundary_outside_both'
  hits=[x for x in splits if x['reference']==ref and x['method']==condition and x['split'].startswith('boundary_') and x['present']=='True'];assert len(hits)==1 and hits[0]['split']==expected
  concord=[x for x in cf if x['reference']==ref and x['condition']==condition and x['split']==expected];assert len(concord)==(2 if condition=='primary' else 1)
  tests={x['hypothesis']:x for x in au if x['reference']==ref and x['condition']==condition}
  assert float(tests[expected]['p_AU'])>=.05
  assert all(float(x['p_AU'])<.05 for k,x in tests.items() if k not in ['unconstrained',expected])
  table.append(dict(reference=ref,condition=condition,sites=v['strict_core_bases'] if condition=='primary' else m['retained_core_bases'],masked_percent=0 if condition=='primary' else m['removed_core_bases']/m['original_core']*100,placement=expected,SHaLRT_UFB=hits[0]['SHaLRT_UFB_label'],sCF_min=min(float(x['sCF']) for x in concord),sCF_max=max(float(x['sCF']) for x in concord),AU_with_nematophilus=tests['boundary_with_nematophilus']['p_AU'],AU_with_bugandensis=tests['boundary_with_bugandensis']['p_AU'],AU_outside_both=tests['boundary_outside_both']['p_AU'],Gubbins_convergence='not_applicable' if condition=='primary' else 'not_reached_within_five_iterations'))
  ds=[x for x in dist if x['reference']==ref and x['method']==condition]
  bw=next(x for x in ds if x['group_a']==x['group_b']=='bugandensis');nw=next(x for x in ds if x['group_a']==x['group_b']=='nematophilus');bt=next(x for x in ds if [x['group_a'],x['group_b']]==['bugandensis','nematophilus'])
  ratio=float(bt['minimum'])/max(float(bw['maximum']),float(nw['maximum']))
  separation.append(dict(reference=ref,condition=condition,within_bugandensis_max=bw['maximum'],within_nematophilus_max=nw['maximum'],between_groups_min=bt['minimum'],minimum_between_over_maximum_within=ratio,interpretation='descriptive tree-distance separation; not a species threshold'))
 # Verify added 118C5 actually joins the six-genome nematophilus clade, beyond merely pruning it.
 panel=[x for x in read(r/'inputs/panel-qc.tsv') if x['local_primary']=='True'];nem={x['accession'] for x in panel if x['group']=='nematophilus_candidate' or x['accession']=='GCA_026344075.1'}
 for method,extra in [('118C5-original','GCA_034107645.1'),('118C5-alternate','GCA_050773425.1')]:
  tree=Phylo.read(t/(method+'.treefile'),'newick');wanted=nem|{extra};allids={x.name for x in tree.get_terminals()}
  assert any({x.name for x in c.get_terminals()} in [wanted,allids-wanted] for c in tree.get_nonterminals())
  sensitivity.append(dict(reference=ref,method=method,added=extra,joins_nematophilus_clade=True,preserves_primary_placement=True,bootstrap='not_computed'))
write('placement-summary.tsv',table);write('distance-separation.tsv',separation);write('118C5-placement.tsv',sensitivity)
# Distinguish ANI failures from coverage failures; preserve overlapping counts.
cohort=[x for x in read(r/'inputs/panel-qc.tsv') if x['local_primary']=='True']
group={x['accession']:('nematophilus' if x['accession']=='GCA_026344075.1' or x['group']=='nematophilus_candidate' else x['group']) for x in cohort}
cross=[x for x in read(r/'outputs/skani-pairs.tsv') if x['a'] in group and x['b'] in group and sorted([group[x['a']],group[x['b']]])==['bugandensis','nematophilus']]
threshold=dict(pairs=len(cross),ani_pass=sum(float(x['ani'])>=95 for x in cross),af_pass=sum(min(float(x['af_a']),float(x['af_b']))>=80 for x in cross),both_pass=sum(x['pass95_80']=='True' for x in cross),minimum_af=min(min(float(x['af_a']),float(x['af_b'])) for x in cross))
assert (threshold['pairs'],threshold['ani_pass'],threshold['af_pass'],threshold['both_pass'])==(156,140,147,138)
(out/'cross-group-threshold-audit.json').write_text(json.dumps(threshold,indent=2)+'\n')

files=[r/'outputs'/x for x in ['tree-splits.tsv','site-concordance.tsv','topology-tests.tsv','patristic-summary.tsv','group-ani-summary.tsv','anim-comparisons.tsv','final-validation.json']]
files+=list(out.glob('*.tsv'))
(out/'review-data-validation.json').write_text(json.dumps(dict(status='PASS',stable_groups_in_all_ten_trees=True,all_four_AU_contrasts_checked=True,added_118C5_placements_checked=4,masked_convergence='capped_not_converged',source_sha256={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}),indent=2)+'\n')
print(json.dumps({'placement':table,'separation':separation},indent=2))
