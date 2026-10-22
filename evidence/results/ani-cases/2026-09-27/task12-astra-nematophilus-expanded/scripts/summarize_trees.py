"""Check trees, quantify named splits and distances, and collect topology tests."""
from pathlib import Path
import csv,json,itertools,re,statistics,hashlib
from Bio import Phylo
r=Path(__file__).resolve().parents[1];out=r/'outputs'
def read(p):return list(csv.DictReader((l for l in open(p) if not l.startswith('#')),delimiter='\t'))
def write(p,rs):
 with open(p,'w') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
rows=[x for x in read(r/'inputs/panel-qc.tsv') if x['local_primary']=='True'];ids={x['accession'] for x in rows};group={x['accession']:('nematophilus' if x['accession']=='GCA_026344075.1' or x['group']=='nematophilus_candidate' else x['group']) for x in rows};groups={g:{q for q in ids if group[q]==g} for g in set(group.values())}
assert {g:len(v) for g,v in groups.items()}=={'bugandensis':26,'nematophilus':6,'boundary':2,'context':4}
def canonical(s,allids):return tuple(sorted(s if len(s)<len(allids-s) else allids-s if len(s)>len(allids-s) else min(tuple(sorted(s)),tuple(sorted(allids-s)))))
def splits(tree):
 allids={x.name for x in tree.get_terminals()};result={}
 for c in tree.find_clades(order='postorder'):
  s={x.name for x in c.get_terminals()}
  if 1<len(s)<len(allids)-1:result[canonical(s,allids)]=c
 return result
wanted={**groups,'boundary_with_nematophilus':groups['boundary']|groups['nematophilus'],'boundary_with_bugandensis':groups['boundary']|groups['bugandensis'],'boundary_outside_both':groups['nematophilus']|groups['bugandensis']}
summary=[];distances=[];cfrows=[];au=[];valid=[]
for reference in ['GCA_019046905.1','GCA_026344075.1']:
 o=r/'local'/reference;t=o/'tests';assert json.loads((t/'validation.json').read_text())['status']=='PASS'
 paths={'primary':o/'tree/primary.treefile','masked':t/'masked.treefile','RY':t/'RY.treefile','118C5-original':t/'118C5-original.treefile','118C5-alternate':t/'118C5-alternate.treefile'}
 for method,path in paths.items():
  tree=Phylo.read(path,'newick');actual={x.name for x in tree.get_terminals()};expected=ids|({'GCA_034107645.1'} if method=='118C5-original' else {'GCA_050773425.1'} if method=='118C5-alternate' else set());assert actual==expected
  assert all(c.branch_length is None or c.branch_length>=0 for c in tree.find_clades())
  # After pruning the single added diagnostic assembly, report the same 38-genome splits.
  if actual!=ids:
   for q in actual-ids:tree.prune(q)
  sp=splits(tree)
  for name,ss in wanted.items():
   c=sp.get(canonical(ss,ids));label=(c.name or str(c.confidence or '')) if c else ''
   summary.append(dict(reference=reference,method=method,split=name,present=c is not None,branch_length=c.branch_length if c else '',SHaLRT_UFB_label=label if method in ['primary','masked','RY'] else 'not_computed',tips=len(actual)))
  if method in ['primary','masked']:
   for a,b in itertools.combinations(sorted(ids),2):distances.append(dict(reference=reference,method=method,a=a,b=b,group_a=group[a],group_b=group[b],patristic_distance=tree.distance(a,b)))
  valid.append(dict(path=str(path),tips=len(actual),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 for condition,seeds in [('primary',[20260927,20260928]),('masked',[20260927])]:
  for seed in seeds:
   prefix=t/(f'primary.scf-{seed}' if condition=='primary' else 'masked.scf')
   tree=Phylo.read(str(prefix)+'.cf.branch','newick');sp=splits(tree);data={str(x['ID']):x for x in read(str(prefix)+'.cf.stat')}
   assert 'BETTER TREE FOUND' not in Path(str(prefix)+'.log').read_text()
   base=Phylo.read(paths[condition],'newick');assert set(sp)==set(splits(base))
   for name,ss in wanted.items():
    c=sp.get(canonical(ss,ids))
    if c is None:continue
    key=c.name if c.name is not None else str(int(c.confidence));x=data[key]
    cfrows.append(dict(reference=reference,condition=condition,seed=seed,split=name,branch_id=key,sCF=x['sCF'],sDF1=x['sDF1'],sDF2=x['sDF2'],informative_sites_per_quartet=x['sN']))
 for condition in ['primary','masked']:
  folder=t/('topology-'+condition)
  if condition=='masked' and (t/'topology-masked-recovery/command.json').exists():
   assert json.loads((t/'topology-masked-recovery/command.json').read_text())['exit_code']==0
   folder=t/'topology-masked-recovery'
  text=(folder/'AU.iqtree').read_text();order=json.loads((folder/'candidate-order.json').read_text());parsed=[]
  # IQ-TREE table has signs after each probability; preserve full row and named AU field.
  lines=text.splitlines();header=next(i for i,l in enumerate(lines) if 'Tree' in l and 'logL' in l and 'p-AU' in l)
  for l in lines[header+1:]:
   tok=l.split()
   if not tok or not tok[0].isdigit():
    if parsed:break
    continue
   clean=[v for v in tok if v not in ['+','-']];num=int(clean[0]);assert 1<=num<=4
   parsed.append(dict(reference=reference,condition=condition,index=num,hypothesis=order[num-1]['hypothesis'],logL=float(clean[1]),deltaL=float(clean[2]),p_AU=float(clean[-1]),raw_row=l.strip()))
  assert len(parsed)==4,(condition,parsed);au+=parsed
write(out/'tree-splits.tsv',summary);write(out/'patristic-pairs.tsv',distances);write(out/'site-concordance.tsv',cfrows);write(out/'topology-tests.tsv',au)
stat=[]
for ref,method in itertools.product(['GCA_019046905.1','GCA_026344075.1'],['primary','masked']):
 for a,b in itertools.combinations_with_replacement(sorted(groups),2):
  vals=[x['patristic_distance'] for x in distances if x['reference']==ref and x['method']==method and sorted([x['group_a'],x['group_b']])==[a,b]]
  stat.append(dict(reference=ref,method=method,group_a=a,group_b=b,pairs=len(vals),minimum=min(vals),median=statistics.median(vals),maximum=max(vals)))
write(out/'patristic-summary.tsv',stat)
(out/'tree-validation.json').write_text(json.dumps({'status':'PASS','trees':valid,'split_rows':len(summary),'site_concordance_rows':len(cfrows),'AU_rows':len(au),'pair_distances':len(distances),'fixed_scf_topologies_verified':True,'no_rank_decision':True},indent=2)+'\n')
print(json.dumps({'AU':au,'boundary_splits':[x for x in summary if x['split'].startswith('boundary_') and x['present']]},indent=2))
