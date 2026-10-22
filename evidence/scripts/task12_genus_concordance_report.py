#!/usr/bin/env python3
"""Join sCF to focal IQ-TREE branches by exact unrooted splits, never row order."""
import argparse,csv,json
from pathlib import Path
import task12_genus_tree_report as m
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run
meta={x['tip_id']:x for x in m.readtsv(r/'inputs/panel.tsv')};rows=[]
for label in ['primary-full80','canonical-masked']:
 t=m.Phylo.read(r/f'trees/{label}-iqtree.treefile','newick');t,_=m.root_tree(t);ids=m.leaves(t.root)
 b=m.Phylo.read(r/f'trees/{label}-scf.cf.branch','newick');assert m.leaves(b.root)==ids
 stats={x['ID']:x for x in csv.DictReader((line for line in open(r/f'trees/{label}-scf.cf.stat') if not line.startswith('#')),delimiter='\t')}
 matches={}
 for n in b.get_nonterminals():
  if n is b.root:continue
  key=m.splitkey(m.leaves(n),ids);bid=n.name or str(int(n.confidence));assert key not in matches;matches[key]=stats[bid]
 clades={name:t.common_ancestor(sorted(pair)) for name,pair in m.PAIRS.items()}
 clades['named_Enterobacter_type_MRCA']=t.common_ancestor(m.read_named_types(meta,ids))
 clades['RHBSTW_parent']=t.get_path('GCA_013782625.1')[-2]
 for name,node in clades.items():
  key=m.splitkey(m.leaves(node),ids);stat=matches[key];assert stat['Label']==node.source_support
  rows.append({'alignment':label,'test':name,'descendant_count':len(m.leaves(node)),'SH_aLRT_UFBoot':node.source_support,'branch_id':stat['ID'],**{k:stat[k] for k in ['sCF','sDF1','sDF2','sN']},'descendants':';'.join(sorted(m.leaves(node)))})
m.write(r/'tables/final-focal-concordance.tsv',rows)
print(json.dumps([{k:v for k,v in x.items() if k!='descendants'} for x in rows],indent=2))
