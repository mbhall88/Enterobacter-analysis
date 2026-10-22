#!/usr/bin/env python3
"""Plot measured ANI/coverage changes; retain the figure source table alongside."""
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/enterobacter-focused-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import task12_focused_scope as s
R=s.RUN
labels={r['selected_accession_version']:r['taxon_name'].replace('Enterobacter ','').replace('hormaechei subsp. ','') for r in s.rows(s.ROOT/'manifests/types.tsv') if r['selected_accession_version']}
labels.update({'GCA_019047005.1':'alternate ludwigii type','GCA_013743295.1':'oharae neighbour'})
rows=s.rows(R/'anim-sensitivity.tsv');(R/'figures').mkdir(exist_ok=True)
fig,axs=plt.subplots(2,3,figsize=(13,8),layout='constrained',sharex='col')
for k,(a,title) in enumerate([('GCA_001750725.1','ludwigii type: add 95.3 kb'),('GCA_013731155.1','oharae RHBSTW-00844: add 361.6 kb')]):
 rr=[r for r in rows if r['pair'].startswith(a)];ys=list(range(len(rr)))
 for j,(metric,lab,threshold) in enumerate([('ani','ANIm (%)',95),('target_af','Focal genome AF (%)',80),('comparator_af','Comparator AF (%)',80)]):
  ax=axs[k,j]
  for y,r in zip(ys,rr):
   x,z=float(r['original_'+metric]),float(r['plus_'+metric]);ax.plot([x,z],[y,y],c='0.55',lw=2);ax.scatter(x,y,c='#3b5b92',marker='o',s=32);ax.scatter(z,y,c='#b64f1d',marker='D',s=28)
  ax.axvline(threshold,c='0.5',ls='--',lw=1);ax.set_xlabel(lab);ax.set_yticks(ys,[labels[r['pair'].split('__')[1]] for r in rr] if j==0 else []);ax.invert_yaxis();ax.grid(axis='x',alpha=.15);ax.spines[['top','right']].set_visible(False)
  if j==0:ax.set_title(title,loc='left',fontsize=12)
fig.legend(handles=[Line2D([],[],color='#3b5b92',marker='o',ls='',label='Original scope'),Line2D([],[],color='#b64f1d',marker='D',ls='',label='Add flagged record')],loc='outside upper right',ncol=2)
fig.text(.01,-.025,'Minimum ANIm and minimum AF per genome across both directions. Historic type-anchor labels identify comparators; they do not assert current rank.\nComparator scope is held fixed, including the documented hormaechei reference limitation. Source: anim-sensitivity.tsv.',fontsize=8)
for ext in ['png','pdf']:fig.savefig(R/'figures'/f'ani-scope-sensitivity.{ext}',dpi=180,bbox_inches='tight')
plt.close(fig)

# Render each paired 51-tip tree when final validation exists; omit pending trees.
if (R/'tree-validation.json').exists():
 from Bio import Phylo
 names={r['accession_version']:r['formal_names'].split(';')[0].replace('Enterobacter ','E. ') for r in s.rows(s.ROOT/'manifests/full.tsv')}
 oharae={r['accession'] for r in s.rows(R/'tree-panel.tsv') if r['oharae_member']=='True'}
 placements={(r['alignment'],r['scope']):r for r in s.rows(R/'tree-placement.tsv')}
 for mask in ['primary','canonical']:
  fig,axes=plt.subplots(1,2,figsize=(21,14),layout='constrained')
  for ax,mode in zip(axes,['original','plus']):
   t=Phylo.read(R/'trees'/f'{mask}-{mode}'/'tree.treefile','newick');t.root_with_outgroup('GCA_001750725.1');t.ladderize()
   def lab(c):
    if not c.is_terminal():return None
    a=c.name
    if a=='GCA_013731155.1':return 'RHBSTW-00844 / '+a+' [FOCAL]'
    return a+' / '+('oharae' if a in oharae else names.get(a,'') or labels.get(a,''))
   colors={lab(c):('#b64f1d' if c.name=='GCA_013731155.1' else '#3265a0' if c.name in oharae else 'black') for c in t.get_terminals()}
   plt.rcParams['font.size']=8
   Phylo.draw(t,axes=ax,label_func=lab,label_colors=colors,do_show=False,show_confidence=False)
   m=placements[(mask,mode)];state='18-member oharae clade' if m['oharae_monophyletic']=='True' else 'oharae not separately resolved'
   ax.set_title(f'{mask.capitalize()} / {mode}\n{state}; encompassing node {m["oharae_support"]}',fontsize=12)
   xmax=max(t.depths().values());ax.set_xlim(-.001,xmax*1.9);ax.set_ylabel('');ax.set_xlabel('Substitutions per amino-acid site');ax.spines[['top','right','left']].set_visible(False);ax.set_yticks([])
  fig.suptitle('Focused scope sensitivity: original chromosome versus chromosome + CP056156.1',fontsize=14)
  fig.text(.015,-.045,'Blue: previously reviewed oharae subgroup. Orange: focal assembly. Node support: SH-aLRT / ultrafast bootstrap. Same51 tips and fixed alignment columns.\nRooted on ludwigii for display only. Historic anchor names identify genomes, not a new rank decision. Identical masked rows cannot resolve some subgroups.',fontsize=10)
  for ext in ['png','pdf']:fig.savefig(R/'figures'/f'tree-{mask}-paired.{ext}',dpi=160,bbox_inches='tight')
  plt.close(fig)
