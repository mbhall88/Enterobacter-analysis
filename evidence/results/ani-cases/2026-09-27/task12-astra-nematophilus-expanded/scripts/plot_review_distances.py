"""Plot quantitative within/between lineage distances; no species cutoff implied."""
from pathlib import Path
import csv,json,statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((r/'outputs/patristic-pairs.tsv').open(),delimiter='\t'))
refs=['GCA_019046905.1','GCA_026344075.1'];rng=np.random.default_rng(20260927)
fig,axes=plt.subplots(2,2,figsize=(12,7.5),sharex=True,sharey=True)
groups=[('bugandensis','bugandensis'),('nematophilus','nematophilus'),('bugandensis','nematophilus')]
labels=['Within bugandensis (325 pairs)','Within nematophilus (15 pairs)','Between groups (156 pairs)'];colours=['#0072B2','#D55E00','#444444'];source=[]
for i,ref in enumerate(refs):
 for j,method in enumerate(['primary','masked']):
  ax=axes[i,j];stats=[]
  for y,((a,b),label,col) in enumerate(zip(groups,labels,colours)):
   subset=[x for x in rows if x['reference']==ref and x['method']==method and sorted([x['group_a'],x['group_b']])==[a,b]];v=[float(x['patristic_distance']) for x in subset];assert len(v)==[325,15,156][y]
   source.extend(subset);ax.scatter(v,y+rng.uniform(-.12,.12,len(v)),s=10,alpha=.25,color=col)
   ax.plot([min(v),max(v)],[y,y],lw=1.2,color=col);ax.plot([statistics.median(v)]*2,[y-.17,y+.17],lw=2.3,color=col);stats.append(v)
  within=max(max(stats[0]),max(stats[1]));between=min(stats[2]);assert between>within
  ax.axvspan(within,between,color='#e8eddc',alpha=.6)
  ax.set_title(('bugandensis' if i==0 else 'nematophilus')+' reference · '+('unmasked' if j==0 else 'masked*'),fontsize=11)
  ax.text(.98,.02,f'Minimum between / maximum within = {between/within:.2f}',transform=ax.transAxes,ha='right',fontsize=8)
  ax.set_yticks(range(3),labels);ax.invert_yaxis();ax.set_xlim(0,.076);ax.grid(axis='x',alpha=.15)
  if i==1:ax.set_xlabel('Patristic distance (sum of branch lengths; substitutions/site)')
  for side in ['top','right']:ax.spines[side].set_visible(False)
fig.suptitle('Quantifying lineage separation in the sampled genomes',fontsize=16)
fig.text(.5,.025,'Dots are genome pairs, not independent replicates. Bars: observed range; ticks: median. Shading: within/between gap.\nRatios describe separation in these trees, not a species threshold. *Masking reached the five-iteration cap without convergence.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.085,1,.94])
for ext in ['png','pdf']:fig.savefig(r/'figures'/('tree-distance-separation.'+ext),dpi=180)
with (r/'review/tree-distance-figure-source.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(source[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(source)
(r/'review/tree-distance-figure-validation.json').write_text(json.dumps({'status':'PASS','plotted_pairs':len(source),'panels':4,'groups_per_panel':3,'expected_pairs_per_group':[325,15,156],'seed':20260927,'species_threshold':None},indent=2)+'\n')
