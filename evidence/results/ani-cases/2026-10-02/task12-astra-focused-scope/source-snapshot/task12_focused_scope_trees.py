#!/usr/bin/env python3
"""Read matched IQ-TREE outputs, test focal clades and draw paired context views."""
import json,os
from pathlib import Path
from Bio import Phylo
import task12_focused_scope as s
R=s.RUN

def splits(tree,positive=False):
 allnames=frozenset(t.name for t in tree.get_terminals());out=set()
 for c in tree.get_nonterminals():
  a=frozenset(t.name for t in c.get_terminals());b=allnames-a
  if min(len(a),len(b))<2 or (positive and (c.branch_length or 0)<1e-8):continue
  out.add(tuple(sorted(min((a,b),key=lambda q:(len(q),tuple(sorted(q)))))))
 return out

def main():
 panel=s.rows(R/'tree-panel.tsv');ids={x['accession'] for x in panel};oharae={x['accession'] for x in panel if x['oharae_member']=='True'};target='GCA_013731155.1'
 output=[];trees={};comparisons=[]
 for mask in ['primary','canonical']:
  for mode in ['original','plus']:
   key=f'{mask}-{mode}';p=R/'trees'/key/'tree.treefile';t=Phylo.read(p,'newick');assert {x.name for x in t.get_terminals()}==ids
   trees[key]=t;t.root_with_outgroup('GCA_001750725.1')
   c=t.common_ancestor(sorted(oharae));members={x.name for x in c.get_terminals()};assert target in members
   leaf=next(x for x in t.get_terminals() if x.name==target);parent=next(c for c in t.find_clades() if leaf in c.clades)
   sister={x.name for x in parent.get_terminals()}-{target}
   output.append(dict(alignment=mask,scope=mode,tips=len(ids),oharae_monophyletic=members==oharae,oharae_mrca_size=len(members),oharae_support=c.name or c.confidence,focal_parent_size=len(sister)+1,focal_parent_support=parent.name or parent.confidence,focal_sisters=';'.join(sorted(sister)),focal_terminal_branch=leaf.branch_length,tree_sha256=s.sha(p)))
  a,b=trees[mask+'-original'],trees[mask+'-plus'];sa,sb=splits(a),splits(b);pa,pb=splits(a,True),splits(b,True)
  comparisons.append(dict(alignment=mask,unrooted_symmetric_split_difference=len(sa^sb),positive_branch_symmetric_split_difference=len(pa^pb),original_nontrivial_splits=len(sa),plus_nontrivial_splits=len(sb)))
 s.write(R/'tree-placement.tsv',output);s.write(R/'tree-split-comparison.tsv',comparisons)
 (R/'tree-validation.json').write_text(json.dumps(dict(trees=4,all_tip_sets_correct=True,focal_membership_tests=output,split_comparisons=comparisons),indent=2))
 print(json.dumps(output,indent=2))

if __name__=='__main__':main()
