import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,TwoSlopeNorm
r=Path(__file__).resolve().parents[1];figdir=r/'figures';figdir.mkdir(exist_ok=True)
read=lambda p:list(csv.DictReader(open(p),delimiter='\t'))
rows=[x for x in read(r/'inputs/panel-qc.tsv') if x['local_primary']=='True'];g={x['accession']:('nematophilus_candidate' if x['accession']=='GCA_026344075.1' else x['group']) for x in rows}
order=['bugandensis','nematophilus_candidate','boundary','context'];ids=sorted(g,key=lambda x:(order.index(g[x]),x));ix={q:i for i,q in enumerate(ids)}
ps=[x for x in read(r/'outputs/skani-pairs.tsv') if x['a'] in g and x['b'] in g];assert len(ps)==703
matrix=np.eye(len(ids))*100
for x in ps:matrix[ix[x['a']],ix[x['b']]]=matrix[ix[x['b']],ix[x['a']]]=float(x['ani'])
colours={'bugandensis':'#3574a8','nematophilus_candidate':'#d47713','boundary':'#982f87','context':'#666666'}
labels={'bugandensis':'bugandensis','nematophilus_candidate':'nematophilus-like','boundary':'ECC-219 + strain 93','context':'chuandaensis context'}
fig,(ax,bx)=plt.subplots(1,2,figsize=(15,7),gridspec_kw={'width_ratios':[1.1,1]},layout='constrained')
cmap=ListedColormap(np.vstack([plt.cm.Blues_r(np.linspace(.2,.8,128)),plt.cm.YlOrRd(np.linspace(.2,.95,128))]));norm=TwoSlopeNorm(vmin=93,vcenter=95,vmax=100)
im=ax.imshow(matrix,cmap=cmap,norm=norm,interpolation='nearest');centres=[];start=0
for k in order:
 count=sum(g[q]==k for q in ids);centres.append(start+(count-1)/2);start+=count
 if start<len(ids):ax.axhline(start-.5,color='black',lw=.7);ax.axvline(start-.5,color='black',lw=.7)
ax.set_xticks(centres,[labels[k]+'\n(n='+str(sum(g[q]==k for q in ids))+')' for k in order],rotation=35,ha='right',fontsize=9)
ax.set_yticks(centres,[labels[k] for k in order],fontsize=9);ax.set_title('All 703 pairs in the strict 38-genome panel')
cbar=fig.colorbar(im,ax=ax,shrink=.8,ticks=[93,94,95,96,98,100]);cbar.set_label('skani ANI (%) — blue <95; orange/red ≥95');cbar.ax.axhline(95,color='black',lw=2)
classes=[('bugandensis','bugandensis'),('nematophilus_candidate','nematophilus_candidate'),('bugandensis','nematophilus_candidate'),('bugandensis','boundary'),('nematophilus_candidate','boundary'),('boundary','boundary')]
rng=np.random.default_rng(20260927);source=[]
for j,(a,b) in enumerate(classes):
 selected=[x for x in ps if sorted([g[x['a']],g[x['b']]])==sorted([a,b])];values=np.array([float(x['ani']) for x in selected]);assert len(values)
 bx.scatter(values,j+rng.uniform(-.14,.14,len(values)),s=13,alpha=.55,c=('#555555' if a!=b else colours[a]),edgecolors='none')
 bx.plot([values.min(),values.max()],[j,j],color='black',lw=1);bx.plot(np.median(values),j,'|',color='black',ms=15,mew=2)
 for x in selected:source.append(dict(category=labels[a]+' / '+labels[b],**x))
 bx.text(100.2,j,f'n={len(values)}',va='center',fontsize=9)
bx.axvspan(93.7,95,color='#d8e7f5',zorder=-2);bx.axvline(95,color='black',ls='--',lw=1.3)
bx.set_yticks(range(len(classes)),[('within '+labels[a] if a==b else labels[a]+' × '+labels[b]) for a,b in classes],fontsize=9)
bx.invert_yaxis();bx.set_xlim(93.7,100.75);bx.set_xticks([94,95,96,97,98,99,100]);bx.set_xlabel('skani ANI (%)');bx.set_title('Within- and between-group separation\nEach dot is a pair; pairs are not independent replicates')
fig.suptitle('Expanded bugandensis / nematophilus sampling',fontsize=16)
fig.get_layout_engine().set(rect=(0,.07,1,.96));fig.text(.01,.015,'Candidate groups are provisional. 118C5 assemblies and the poor-quality NCTC 13925 record are excluded from this strict panel.\nANI is the lower reciprocal estimate; aligned fractions and method-specific ANIm checks are in the source tables.',fontsize=9)
for ext in ['png','pdf']:fig.savefig(figdir/('expanded-ani.'+ext),dpi=180)
with (figdir/'expanded-ani-source.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(source[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(source)
print('38 genomes, 703 pairs plotted')
