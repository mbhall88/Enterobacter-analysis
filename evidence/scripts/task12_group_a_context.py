#!/usr/bin/env python3
"""Reuse saved trees and branch concordance for population evidence; no tree inference."""
from pathlib import Path
import argparse,csv,json
import numpy as np
from Bio import Phylo
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from task12_group_a_populations import rows,write,BASE,CFG,TYPES,sha
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output
members=rows(out/'outputs/assignments.tsv');groups={r['accession']:int(r['operational_cluster']) for r in members};ids=set(groups)
# Match an unrooted split, never an arbitrary Newick root or the order of branch rows.
records=[];sources=[]
for variant in ['gcf','gcfc']:
 b=BASE/f'core-genome/concordance-gcf-all/ml.{variant}.cf.branch';s=BASE/f'core-genome/concordance-gcf-all/ml.{variant}.cf.stat';sources.extend([b,s]);t=Phylo.read(b,'newick');universe={x.name for x in t.get_terminals()}
 with s.open() as f: stats={r['ID']:r for r in csv.DictReader((l for l in f if not l.startswith('#')),delimiter='\t')}
 for c in sorted(set(groups.values())):
  target={q for q in universe & ids if groups[q]==c};found=[]
  for node in t.get_nonterminals():
   ds={x.name for x in node.get_terminals()}
   if ds==target or universe-ds==target:
    bid=str(int(node.confidence)) if node.confidence is not None else node.name
    if bid in stats:found.append(stats[bid])
  assert len(found)<=1
  r=found[0] if found else {}
  records.append(dict(cluster=c,represented=len(target),gene_tree_variant='raw' if variant=='gcf' else 'collapsed_below_one_expected_substitution',reference_split_present=bool(found),branch_id=r.get('ID',''),gene_concordance_percent=r.get('gCF',''),concordant_gene_trees=r.get('gCF_N',''),decisive_gene_trees=r.get('gN','')))
write(out/'outputs/existing-gene-concordance.tsv',records)
colors={1:'#0072B2',2:'#777777',3:'#E69F00',4:'#009E73',5:'#CC79A7',6:'#56B4E9',7:'#D55E00'}
names={1:'nototheniae-type population',2:'mixed hoffmannii-type cluster',3:'population P3',4:'intestinihominis-type population',5:'population P5',6:'population P6',7:'population P7'}
cfg=json.loads(CFG.read_text());fig,axes=plt.subplots(1,3,figsize=(20,13))
for ax,(method,v) in zip(axes,cfg['comparison_sets']['ingroup-1076']['rooting']['per_tree'].items()):
 path=Path(v['rootings']['both_set12_outgroups']['tree']);sources.append(path);tree=Phylo.read(path,'newick');root=tree.common_ancestor(sorted(ids));assert {x.name for x in root.get_terminals()}==ids
 # Collapse only exact, monophyletic operational clusters; preserve residual C2 genomes individually.
 terminals=[];kids={};desc={}
 def walk(node):
  ss={x.name for x in node.get_terminals()};desc[node]=ss;cs={groups[q] for q in ss}
  collapse=len(cs)==1 and len(ss)>1 and ss=={q for q in ids if groups[q]==next(iter(cs))}
  if collapse or node.is_terminal():terminals.append(node);kids[node]=[]
  else:
   kids[node]=node.clades
   for k in node.clades:walk(k)
 walk(root);y={node:len(terminals)-i-1 for i,node in enumerate(terminals)};x={root:0.}
 def calc(node):
  for k in kids[node]:x[k]=x[node]+(k.branch_length or 0);calc(k)
  if kids[node]:y[node]=np.mean([y[k] for k in kids[node]])
 calc(root);maxx=max(x.values());offset=maxx*.025
 for node in kids:
  if kids[node]:
   ax.plot([x[node],x[node]],[min(y[k] for k in kids[node]),max(y[k] for k in kids[node])],color='#444',lw=.8)
   for k in kids[node]:
    cs={groups[q] for q in desc[k]};color=colors[next(iter(cs))] if len(cs)==1 else '#444';ax.plot([x[node],x[k]],[y[k],y[k]],color=color,lw=1.3)
 for node in terminals:
  ss=desc[node];c=groups[next(iter(ss))];types=[TYPES[q] for q in ss if q in TYPES]
  if len(ss)>1:label=f'P{c}: {names[c]} (n={len(ss)})';ax.scatter(x[node],y[node],marker='>',s=65,color=colors[c])
  else:
   q=next(iter(ss));label=TYPES.get(q,q)+' [P'+str(c)+']';ax.scatter(x[node],y[node],s=15,color=colors[c])
  ax.text(x[node]+offset,y[node],label,va='center',fontsize=7.5,fontweight='bold' if types else 'normal',color=colors[c])
 ax.set_ylim(-1,len(terminals));ax.set_xlim(-offset,maxx*2.5);ax.set_yticks([]);ax.set_xlabel('Original branch length');ax.set_title({'bac120_fasttree':'BAC120 / FastTree','verticall_fastme':'Verticall / FastME','mashtree_consensus':'Mashtree consensus'}[method]);ax.spines[['top','right','left']].set_visible(False)
fig.suptitle('Group A: independent core-population clusters on the three existing global trees\nSix clusters form clades; the grey hoffmannii-type cluster is dispersed. Colours are not subspecies assignments.',fontsize=13)
fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(out/'figures/populations-on-three-trees.png',dpi=170);fig.savefig(out/'figures/populations-on-three-trees.pdf');plt.close(fig)
(out/'context-inputs.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in sources))
print(json.dumps(records,indent=2))
