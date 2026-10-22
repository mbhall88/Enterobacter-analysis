#!/usr/bin/env python3
"""Exact-clade summaries and full/collapsed displays of genus-context trees."""
import argparse,copy,csv,hashlib,json,math
from pathlib import Path
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
OUT={'GCF_001598475.1','GCF_017161565.1'}
FOCAL={'GCA_013782625.1':'RHBSTW-00994','GCA_013898335.1':'RHB15-C17','GCA_014169315.1':'cloacae_P query 1','GCA_015708675.1':'cloacae_P query 2','GCA_043968675.1':'ECC-175'}
PAIRS={'cloacae_P_pair':{'GCA_014169315.1','GCA_015708675.1'},'nimipressuralis_control':{'GCA_013898335.1','GCF_004115925.1'},'ECC175_new_neighbour':{'GCA_043968675.1','GCF_041399425.1'}}
COL={'Enterobacter':'#0072B2','Lelliottia':'#D55E00','Leclercia':'#009E73','Huaxiibacter':'#CC79A7','Silvania':'#8B6508','Pseudenterobacter':'#8064A2'}

def readtsv(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
 if not rows:return
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def leaves(n):return frozenset(x.name for x in n.get_terminals())
def splitkey(side,allids):
 other=allids-side;return tuple(sorted(side if (len(side),tuple(sorted(side)))<(len(other),tuple(sorted(other))) else other))
def support(n):return n.name if n.name and not n.is_terminal() else ('' if n.confidence is None else str(n.confidence))
def genus(row):
 s=row['NCBI_name'].removeprefix('uncultured ').split()[0]
 return s if s in COL else 'context/unassigned'
def root_tree(tree):
 ids=leaves(tree.root);sm={splitkey(leaves(n),ids):support(n) for n in tree.get_nonterminals() if n is not tree.root}
 # Find the exact unrooted outgroup edge, not merely a convenient MRCA.
 candidates=[n for n in tree.find_clades() if leaves(n)==OUT]
 assert len(candidates)==1,'Prespecified outgroup edge not recovered'
 before={(a,b):tree.distance(a,b) for a in sorted(ids) for b in sorted(ids) if a<b}
 tree.root_with_outgroup(candidates[0]);tree.ladderize()
 for (a,b),d in before.items():assert math.isclose(tree.distance(a,b),d,abs_tol=1e-10)
 for n in tree.find_clades():n.source_support=sm.get(splitkey(leaves(n),ids),'')
 return tree,len(before)
def save_rooted(tree,path):
 # Re-rooting moves nodes: serialize the original split support, not the
 # confidence that Bio.Phylo left attached to a formerly different edge.
 export=copy.deepcopy(tree)
 for n in export.get_nonterminals():
  n.name=n.source_support or None;n.confidence=None
 Phylo.write(export,path,'newick',format_branch_length='%1.12f')
 restored=Phylo.read(path,'newick');ids=sorted(leaves(tree.root))
 for i,a in enumerate(ids):
  for b in ids[i+1:]:assert math.isclose(tree.distance(a,b),restored.distance(a,b),abs_tol=1e-9)
def brief_types(meta,ids):return sorted(n for n,x in meta.items() if n in ids and x['type_material'] and x['NCBI_name'].startswith('Enterobacter '))
def read_named_types(meta,ids):return brief_types(meta,ids)

def diagram(t,meta,anchors,path,title,collapse):
 protected=set(FOCAL)|set(anchors)|{'GCF_041399425.1','GCF_001654845.1','GCF_019047785.1'}
 displayed=[];children={};collapsed=set()
 def walk(n):
  ls=leaves(n);gs={genus(meta[x]) for x in ls}
  if collapse and len(ls)>=2 and not(ls & protected) and len(gs)==1 and next(iter(gs))!='context/unassigned':
   collapsed.add(n);displayed.append(n);children[n]=[]
  elif n.is_terminal():displayed.append(n);children[n]=[]
  else:
   children[n]=n.clades
   for c in n.clades:walk(c)
 walk(t.root);ys={n:i for i,n in enumerate(displayed)}
 def pos(n):
  if n not in ys:ys[n]=sum(pos(c) for c in children[n])/len(children[n])
  return ys[n]
 pos(t.root);depth=t.depths();mx=max(depth[n] for n in displayed);font=8.5 if not collapse else 10
 fig,ax=plt.subplots(figsize=(18,max(7,len(displayed)*(.22 if not collapse else .32))))
 sources=[]
 def draw(n):
  kids=children[n]
  if kids:
   ax.plot([depth[n]]*2,[min(ys[c] for c in kids),max(ys[c] for c in kids)],color='#555555',lw=.75)
   for c in kids:
    ax.plot([depth[n],depth[c]],[ys[c],ys[c]],color='#555555',lw=.75);draw(c)
  else:
   ls=leaves(n)
   if n in collapsed:
    g=genus(meta[next(iter(ls))]);lab=f'{g}: {len(ls)} reference tips collapsed';color=COL.get(g,'#555555')
    ax.scatter([depth[n]],[ys[n]],marker='>',s=40,color=color)
   else:
    x=n.name;row=meta[x];label=row['NCBI_name'];strain=json.loads(row['strain']);st=strain.get('strain',strain.get('isolate',''))
    if x in FOCAL:lab=f'{FOCAL[x]}  |  {x}';color='#A00024';ax.scatter([depth[n]],[ys[n]],s=20,color=color,zorder=5)
    else:
     lab=f'{label}'+(f' [{st}]' if st and st not in label else '')+f'  |  {x}';color=COL.get(genus(row),'#555555')
    if x in anchors:lab+='  [GENUS TYPE]' if anchors[x]['status'].startswith('validly') else '  [PROPOSED TYPE]'
    if row['selection']!='primary':lab+='  [LOW COVERAGE]'
   ax.text(depth[n]+mx*.016,ys[n],lab,fontsize=font,va='center',color=color,fontweight='bold' if ls & set(FOCAL) else 'normal')
  if n is not t.root and kids and (leaves(n)&protected) and n.source_support:
   ax.text(depth[n],ys[n]-.18,n.source_support,fontsize=6.5,color='#555555',ha='right',va='bottom')
  sources.append(dict(display_node=hashlib.sha256('|'.join(sorted(leaves(n))).encode()).hexdigest()[:12],descendant_count=len(leaves(n)),descendants=';'.join(sorted(leaves(n))),x=depth[n],y=ys[n],support=n.source_support,collapsed=n in collapsed))
 draw(t.root)
 ax.set_xlim(-mx*.01,mx*2.8);ax.set_ylim(len(displayed),-1.5);ax.set_yticks([]);ax.set_xlabel('Amino-acid substitutions per site',loc='left')
 detail='Red: focal chromosomes. Names: NCBI source labels, not new decisions. Root: sampled Erwinia/Pantoea clade.'
 support_label='FastTree SH-like local support (0–1)' if 'FastTree' in title else 'SH-aLRT / ultrafast bootstrap (%)'
 ax.set_title(title+'\n'+detail+'\nSupport: '+support_label+'. GTDB crosswalk and all accessions: panel.tsv.',loc='left',fontsize=10,pad=18,linespacing=1.5)
 for sp in ['left','right','top']:ax.spines[sp].set_visible(False)
 ax.spines['bottom'].set_bounds(0,mx)
 ticks=[0,mx/2,mx];ax.set_xticks(ticks,[f'{v:.3f}' for v in ticks]);fig.tight_layout()
 fig.savefig(Path(str(path)+'.png'),dpi=150);fig.savefig(Path(str(path)+'.pdf'));plt.close(fig)
 write(Path(str(path)+'.nodes.tsv'),sources)
 assert sum(len(leaves(n)) for n in displayed)==len(t.get_terminals())
 return len(displayed)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--stage',choices=['preliminary','fasttree-complete','final'],required=True);a=p.parse_args();r=a.run
 meta={x['tip_id']:x for x in readtsv(r/'inputs/panel.tsv')};anchors={x['tip_id']:x for x in readtsv(r/'inputs/type-anchors.tsv')}
 targets=[('primary-full80-fasttree.nwk','Fuller BAC120 · FastTree · 40,254 columns'),('canonical-masked-fasttree.nwk','GTDB-masked BAC120 · FastTree · 5,010 columns'),('candidates-masked-fasttree.nwk','GTDB-masked BAC120 · FastTree · all 97 candidates')]
 if a.stage!='preliminary':targets += [('primary-full80-no-composition-outlier-fasttree.nwk','Fuller BAC120 · FastTree · composition-outlier omission')]
 if a.stage=='final':targets +=[(f'{x}-iqtree.treefile',f'{"Fuller" if x.startswith("primary") else "GTDB-masked"} BAC120 · IQ-TREE') for x in ['primary-full80','canonical-masked']]
 dest=r/'figures'/a.stage;dest.mkdir(exist_ok=False)
 quantitative=[];mrca=[];sisters=[];summaries=[]
 for filename,title in targets:
  f=r/'trees'/filename;t=Phylo.read(f,'newick');ids=leaves(t.root);expected=set(meta) if filename.startswith('candidates') else {n for n,x in meta.items() if x['selection']=='primary'}
  if 'no-composition-outlier' in filename:expected.remove('GCA_023425455.1')
  assert ids==expected
  assert all(n.branch_length is None or n.branch_length>=0 for n in t.find_clades());t,nchecks=root_tree(t)
  named=read_named_types(meta,ids);core=t.common_ancestor(named);coreids=leaves(core)
  for label,pair in PAIRS.items():
   node=t.common_ancestor(sorted(pair));quantitative.append(dict(tree=filename,test=label,exact_clade=leaves(node)==pair,tip_count=len(leaves(node)),support=node.source_support,descendants=';'.join(sorted(leaves(node)))))
  for q in FOCAL:
   path=t.get_path(q);parent=path[-2];sib=set(leaves(parent))-{q}
   sisters.append(dict(tree=filename,query=q,parent_tips=len(leaves(parent)),parent_support=parent.source_support,sister_tips=len(sib),sister_members=';'.join(sorted(sib)),inside_named_Enterobacter_type_MRCA=q in coreids,named_type_count=len(named),type_MRCA_tips=len(coreids),type_MRCA_support=core.source_support))
   for anchor,rec in anchors.items():
    n=t.common_ancestor(q,anchor);mrca.append(dict(tree=filename,query=q,anchor_genus=rec['genus'],anchor=anchor,mrca_tip_count=len(leaves(n)),mrca_support=n.source_support,patristic_distance=t.distance(q,anchor)))
  save_rooted(t,dest/(filename+'.rooted.nwk'))
  overview=diagram(t,meta,anchors,dest/(filename+'.overview'),title+(' · IQ-TREE checks pending' if a.stage!='final' else ''),True)
  full=diagram(t,meta,anchors,dest/(filename+'.full'),title+(' · IQ-TREE checks pending' if a.stage!='final' else ''),False)
  summaries.append(dict(tree=filename,tips=len(ids),distance_checks=nchecks,overview_rows=overview,full_rows=full,named_Enterobacter_type_tips=len(named),named_type_MRCA_size=len(coreids)))
 write(r/'tables'/f'{a.stage}-clade-tests.tsv',quantitative);write(r/'tables'/f'{a.stage}-focal-placement.tsv',sisters);write(r/'tables'/f'{a.stage}-genus-anchor-distances.tsv',mrca)
 (r/'tables'/f'{a.stage}-figure-validation.json').write_text(json.dumps(summaries,indent=2)+'\n');print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
