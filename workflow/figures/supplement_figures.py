"""Render saved type projections and ANI measurements; no phylogenetic inference."""
import copy,csv,itertools,json,math,os,textwrap
from pathlib import Path
from prepare import ROOT,RUN,sha,rows
os.environ['MPLCONFIGDIR']=str(RUN/'outputs/mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
OUT=RUN/'outputs/figures';DATA=OUT/'source-data';DATA.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
TYPE=RUN/'outputs/tables/S02-type-display-map.tsv'
ANI=ROOT/'results/ani-screen/2026-09-17/task04-opus-type-anchored/tables/type-versus-type-ani.tsv'
SUMMARY=ROOT/'results/ani-cases/2026-09-28/task12-astra-ab-review/outputs-v2/summary.json'
GENOMES=SUMMARY.with_name('genomes.tsv')
METHODS=['bac120_fasttree','verticall_fastme','mashtree_consensus']
TREE={k:ROOT/f'results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/type-anchor-projection-35/trees/{k}.nwk' for k in METHODS}
mp=rows(TYPE);byid={i+1:r for i,r in enumerate(mp)}
pairs={frozenset((r['type_a_accession_version'],r['type_b_accession_version'])):r for r in rows(ANI)}
assert len(pairs)==528
titles=['BAC120 / FastTree','Verticall / FastME','Mashtree consensus']
panels=[
 ('E. xiangfangensis and E. hormaechei',[7,8,9,11,28,34,35,10,5], 'The oharae, steigerwaltii, xiangfangensis, hoffmannii, intestinihominis and nototheniae types are shown with the hormaechei type and kobei/asburiae context. Under the proposed treatment, the first six types belong to E. xiangfangensis; the hormaechei type remains in E. hormaechei. The three hoffmannii-associated types have an unresolved internal history (S04).'),
 ('E. cloacae and E. dissolvens',[1,24,13,15,33], 'The cloacae and dissolvens types are shown with ludwigii, sichuanensis and chinensis context. The proposed treatment retains two subspecies within E. cloacae. The alternative placement of E. ludwigii remains unresolved.'),
 ('E. asburiae and E. dykesii',[5,22,10,12,13], 'The asburiae and dykesii types are shown with kobei, roggenkampii and ludwigii context. The proposed treatment places E. dykesii as a later heterotypic synonym of E. asburiae; relationships of the intermediate population are described in S06.'),
 ('E. bugandensis and E. nematophilus',[25,29,14,32,12], 'The bugandensis and nematophilus types are shown with chengduensis, rongchengensis and roggenkampii context. The named cores remain separate species despite overlapping ANI evidence. Type-only projections do not show the unassigned boundary genomes.'),
 ('E. chengduensis and E. rongchengensis',[14,32,25,29,6], 'The chengduensis and rongchengensis types are shown with bugandensis, nematophilus and soli context. Their proposed treatment as subspecies retains the conflicting alignment-coverage and published dDDH evidence in S06.'),
 ('E. roggenkampii and neighbouring types',[12,4,30,21,25], 'The roggenkampii, quasiroggenkampii, pseudoroggenkampii and vonholyi types are shown with bugandensis context. The proposed quasiroggenkampii subspecies treatment and unresolved pseudoroggenkampii/vonholyi boundary are distinct interpretations; close tree placement alone does not assign the boundary genomes.'),
 ('E. mori / E. quasimori and E. sichuanensis / E. chinensis',[2,18,15,33,16], 'The mori and quasimori types, and the sichuanensis and chinensis types, are shown with chuandaensis context. E. quasimori is treated as a later heterotypic synonym of E. mori, whereas E. sichuanensis and E. chinensis remain separate species. Their conflicting tree placements are examined in S04.'),
]
checks=[];manifest=[];captions=[]
def save(fig,n,sources):
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figureS{n}.{ext}',dpi=300,facecolor='white',bbox_inches='tight',pad_inches=.08,metadata={'Creator':'Saved-evidence manuscript rendering'} if ext=='pdf' else None)
 plt.close(fig)
 manifest.append({'figure':f'S{n}','sources':[{'path':str(p),'sha256':sha(p)} for p in sources], 'formats':['png','pdf','svg'],'new_biological_inference':False})
for num,(title,ids,interpretation) in enumerate(panels,1):
 chosen=[byid[i]['accession'] for i in ids];labels={byid[i]['accession']:f'T{i:02}' for i in ids}
 fig=plt.figure(figsize=(7.2,7.8));grid=fig.add_gridspec(3,6,height_ratios=[2.1,3.0,1.8],hspace=.45,wspace=.5)
 geometry={}
 for j,(method,title2) in enumerate(zip(METHODS,titles)):
  original=Phylo.read(TREE[method],'newick');tree=copy.deepcopy(original)
  for tip in list(tree.get_terminals()):
   if tip.name not in chosen:tree.prune(tip)
  assert {t.name for t in tree.get_terminals()}==set(chosen)
  for a,b in itertools.combinations(chosen,2):
   checks.append({'name':f'S{num} {method} path {a} {b}','pass':math.isclose(original.distance(a,b),tree.distance(a,b),abs_tol=1e-10)})
  # Every retained node is from the original tree, with its original confidence.
  original_groups={frozenset(t.name for t in c.get_terminals() if t.name in chosen):c.confidence for c in original.find_clades() if len([t for t in c.get_terminals() if t.name in chosen])>=2}
  terms=tree.get_terminals();dep=tree.depths();yy={c:len(terms)-i for i,c in enumerate(terms)}
  def coord(c):
   if c not in yy:yy[c]=sum(coord(x) for x in c.clades)/len(c.clades)
   return yy[c]
  coord(tree.root);ax=fig.add_subplot(grid[0,j*2:j*2+2]);maximum=max(dep.values())
  def draw(c):
   if c.clades:
    ax.plot([dep[c]]*2,[min(yy[x] for x in c.clades),max(yy[x] for x in c.clades)],color='#333333',lw=.8)
    for child in c.clades:
     ax.plot([dep[c],dep[child]],[yy[child]]*2,color='#333333',lw=.8);draw(child)
    if c.confidence is not None and c.confidence>=.95:ax.plot(dep[c],yy[c],'o',color='#333333',ms=3)
   else:ax.text(dep[c]+maximum*.025,yy[c],labels[c.name],fontsize=9,va='center')
  draw(tree.root);ax.set_ylim(0,len(terms)+1);ax.set_xlim(0,maximum*1.35);ax.set_yticks([])
  ax.set_title(f'({chr(97+j)}) {title2}',fontsize=9,loc='left');ax.spines['left'].set_visible(False)
  ax.tick_params(axis='x',labelsize=7);ax.locator_params(axis='x',nbins=3)
  ax.set_xlabel('Source distance units',fontsize=7)
  geometry[method]={'tips':[x.name for x in terms],'branches':[{'name':c.name,'length':c.branch_length,'confidence':c.confidence,'descendants':sorted(t.name for t in c.get_terminals())} for c in tree.find_clades()]}
  Phylo.write(tree,DATA/f'figureS{num}-{method}.nwk','newick',format_branch_length='%1.12g',format_confidence='%1.12g')
 source_rows=[]
 for i,a in enumerate(chosen):
  for b in chosen[i+1:]:source_rows.append({'figure':f'S{num}','display_a':labels[a],'display_b':labels[b],**pairs[frozenset((a,b))]})
 with (DATA/f'figureS{num}-pairs.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(source_rows[0]),delimiter='\t');w.writeheader();w.writerows(source_rows)
 for panel,key,cmap,vmin,vmax,label in [(0,'ani_percent','viridis',85,100,'ANI (%)'),(1,'af_min_percent','cividis',0,100,'Minimum AF (%)')]:
  matrix=np.full((len(ids),len(ids)),np.nan)
  for i,a in enumerate(chosen):
   for j,b in enumerate(chosen):
    if i!=j:matrix[i,j]=float(pairs[frozenset((a,b))][key])
  assert np.nanmin(matrix)>=vmin and np.nanmax(matrix)<=vmax
  ax=fig.add_subplot(grid[1,panel*3:panel*3+3]);cm=plt.get_cmap(cmap).copy();cm.set_bad('#eeeeee')
  im=ax.imshow(matrix,cmap=cm,vmin=vmin,vmax=vmax)
  ax.set_xticks(range(len(ids)),[f'T{i:02}' for i in ids],fontsize=7,rotation=45);ax.set_yticks(range(len(ids)),[f'T{i:02}' for i in ids],fontsize=7)
  for i in range(len(ids)):
   for j in range(len(ids)):
    v=matrix[i,j];s='—' if i==j else f'{v:.1f}'
    norm=(v-vmin)/(vmax-vmin);ax.text(j,i,s,ha='center',va='center',fontsize=6.5 if len(ids)>6 else 7.3,color='white' if not math.isnan(v) and norm<.52 else 'black')
  ax.set_title(f'({chr(100+panel)}) {label}',fontsize=9,loc='left');fig.colorbar(im,ax=ax,fraction=.035,pad=.03).ax.tick_params(labelsize=7)
 keyax=fig.add_subplot(grid[2,:]);keyax.axis('off')
 for i,n in enumerate(ids):
  r=byid[n];name=r['historical_type_names'].replace('Enterobacter ','E. ')
  # Collapse redundant shared type aliases for the display key; the full aliases remain Table S03.
  aliases=name.split(';');name=next((x for x in aliases if 'subsp.' in x),aliases[0])
  keyax.text(0,1-i*.115,f'T{n:02}  {name}  ({r["accession"]})',fontsize=8,va='top')
 fig.suptitle(f'Figure S{num}. {title}',fontsize=11,y=.985)
 fig.subplots_adjust(left=.075,right=.97,top=.90,bottom=.035)
 (DATA/f'figureS{num}-tree-geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
 save(fig,num,[ANI,TYPE,*TREE.values()])
 captions.append(f'### Figure S{num}. {title}\n\n![Figure S{num}](figures/figureS{num}.png)\n\n{interpretation} Panels (a–c) show the three saved type-tree projections; panels (d–e) show skani ANI and the smaller of the two aligned fractions. Labels identify the historical type names, not automatic classifications of surrounding genomes. The common display rules above apply.\n<!-- sources: FIGSUP -->\n')

# Draw saved summary quantiles, not new estimates or artificial independent samples.
summary=json.loads(SUMMARY.read_text());selected=[x for x in summary['summaries'] if x['comparison'].startswith('union:')]
assert [x['pairs'] for x in selected]==[14878,95703,75774,2249,5694]
labels=['Within hoffmannii lineage','Within oharae / steigerwaltii / xiangfangensis set','Between the two E. xiangfangensis sets','hoffmannii lineage vs E. hormaechei','Other E. xiangfangensis set vs E. hormaechei']
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,4.8),gridspec_kw={'width_ratios':[2.3,1]},layout='constrained')
boxes=[{'label':str(i+1),'med':s['median'],'q1':s['q25'],'q3':s['q75'],'whislo':s['min'],'whishi':s['max'],'fliers':[]} for i,s in enumerate(selected)]
a.bxp(boxes,positions=range(5),orientation='horizontal',showfliers=False,patch_artist=True,boxprops={'facecolor':'#56B4E9'},medianprops={'color':'black'});a.set_yticks(range(5),['1','2','3','4','5']);a.invert_yaxis();a.set_xlabel('skani ANI (%)');a.set_xlim(93.5,100.2);a.axvline(95,color='#888888',ls=':',lw=1);a.set_title('(a) Pairwise ANI',loc='left',fontsize=10)
fractions=[100*s['passes_95_80']/s['pairs'] for s in selected]
b.barh(range(5),fractions,color='#0072B2');b.invert_yaxis();b.set_yticks([]);b.set_xlim(0,120);b.set_xticks([0,50,100]);b.set_xlabel('Pairs passing (%)');b.set_title('(b) ANI / AF screen',loc='left',fontsize=10)
for i,(s,v) in enumerate(zip(selected,fractions)):b.text(v+2,i,f'{v:.1f}',va='center',fontsize=8)
fig.suptitle('Figure S8. Population comparisons within the broader hormaechei-related clade',fontsize=10)
save(fig,8,[SUMMARY,GENOMES])
with (DATA/'figureS8-summary.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['display','description',*selected[0]],delimiter='\t');w.writeheader();w.writerows({'display':i,'description':labels[i-1],**s} for i,s in enumerate(selected,1))
key='; '.join(f'{i+1}, {label[0].lower()+label[1:]} (n = {s["pairs"]:,} pairs)' for i,(label,s) in enumerate(zip(labels,selected)))
captions.append(f'### Figure S8. Population ANI and alignment coverage\n\n![Figure S8](figures/figureS8.png)\n\nThe numbered comparisons are: {key}. The sets contain 173 hoffmannii-lineage genomes, 438 oharae/steigerwaltii/xiangfangensis-associated genomes and thirteen E. hormaechei genomes. Boxes show the median and interquartile range; whiskers show the minimum and maximum, not confidence intervals. Panel (b) gives the proportion passing ANI ≥95% and AF ≥80% for both genomes. The 438-genome set includes one peripheral genome assigned only to species. Pairwise observations share genomes and are not independent replicates. These are the saved chromosome-scope skani comparisons, not whole-type or reciprocal ANIm measurements.\n<!-- sources: C11,GROUPB,FIGSUP -->\n')
common='''## Supplementary figures

Figures S1–S7 show selected type strains and named neighbours pruned from the same three saved 35-tip projections used for main Figure 2. They orient the taxonomic comparisons; they do not replace the population analyses or establish that omitted genomes lack conflicting relationships. No tree was re-inferred or rerooted. Pairwise path lengths among retained tips are preserved. Distance scales are method-specific and should not be compared across panels. Dots retain source SH-like local support (BAC120) or sketch-seed frequency (Mashtree) ≥0.95 for the original branches; support was not recalculated after pruning, and FastME has no branch-support estimates.

The heatmaps use the original skani type–type comparisons of selected chromosome inputs, including all contigs for draft types. Minimum AF is the smaller aligned fraction of the two inputs, not a reciprocal ANIm coverage result. Diagonals are unmeasured self comparisons and are left blank. Colour scales are fixed across all seven figures. Cell labels are rounded to one decimal place; exact values are supplied in the figure source data. Historical type names identify the sequences; current proposed names and full aliases are given in Table S03. Figure S8 adds the population-level ANI evidence for E. xiangfangensis and E. hormaechei. Standalone PDF, SVG and 300 dpi PNG figures and their exact source tables accompany the supplement.
<!-- sources: FIGSUP,C04,ANI04 -->

'''
(RUN/'outputs/supplement-figures.md').write_text(common+'\n'.join(captions))
(OUT/'supplement-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(RUN/'checks/supplement-figures.json').write_text(json.dumps({'checks':len(checks),'passed':sum(x['pass'] for x in checks),'failed':[x for x in checks if not x['pass']],'details':checks},indent=2)+'\n')
assert all(x['pass'] for x in checks)
print('Rendered eight supplementary figures;',len(checks),'retained path-length checks passed.')
