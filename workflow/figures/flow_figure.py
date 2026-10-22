"""Three-panel manuscript flowchart; every cohort count is bound to saved manifests."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def draw(save,RUN):
    d=json.loads((RUN/'outputs/dataset-flow.json').read_text());n=d['counts']
    fig=plt.figure(figsize=(7.5,9.5))
    a=fig.add_axes([.025,.665,.95,.31]);b=fig.add_axes([.025,.16,.455,.47]);c=fig.add_axes([.52,.16,.455,.47])
    for ax in (a,b,c):ax.set_xlim(-1,101);ax.set_ylim(0,100);ax.axis('off')
    def box(ax,x,y,w,h,s,color='#F1F3F5',fs=8.8):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.25,rounding_size=1',facecolor=color,edgecolor='#444444',lw=.8))
        ax.text(x+w/2,y+h/2,s,ha='center',va='center',fontsize=fs,linespacing=1.25)
    def arrow(ax,p,q):ax.annotate('',xy=q,xytext=p,arrowprops={'arrowstyle':'-|>','lw':.9,'color':'#333333'})
    a.text(0,100,'(a)  Study collection and phylogenetic datasets',fontweight='bold',fontsize=10,va='top')
    box(a,0,64,28,23,f"Initial NCBI collection\nn = {n['initial_collection']:,}")
    box(a,36,64,28,23,f"Study collection\nn = {n['study']:,}", '#DCEEF7')
    box(a,72,64,28,23,f"Classification comparison\nn = {n['classification']:,}", '#DCEEF7',8.4)
    arrow(a,(28,75),(36,75));arrow(a,(64,75),(72,75))
    a.text(32,60,'5 removed\n4 duplicate submissions\n1 unsuitable assembly',ha='center',va='top',fontsize=7.5)
    a.text(83,60,'4 later sequence-quality\nexclusions (Supplement S01)',ha='center',va='top',fontsize=7.5)
    box(a,0,1,100,28,'Phylogenetic inputs: 1,072 study + 16 additional types + 2 outgroups = 1,090\n'
        'Two study quality exclusions preceded inference; a third was pruned → 1,089\n'
        'Ingroup comparison: 1,076 tips; later quality-exclusion sensitivity: 1,088 tips',fs=8.4)
    arrow(a,(50,64),(50,31))
    b.text(0,100,'(b)  References and study evaluation',fontweight='bold',fontsize=10,va='top')
    box(b,0,79,100,14,f"Classifier source set: n = {n['reference_selection']:,}\n1,070 study + 16 additional types", '#DCEEF7')
    box(b,0,53,47,18,f"Primary references\nn = {n['primary_references']}\n124 study + 16 types",fs=8.6)
    box(b,53,53,47,18,f"Study queries\nn = {n['study_queries']}\n922 named + 24 other",fs=8.6)
    arrow(b,(25,79),(25,71));arrow(b,(77,79),(77,71))
    box(b,0,28,100,17,'Separate genomovar component: n = 173\n172 study genomes + 1 additional type\nDefinition members overlap the source set',fs=8.6)
    box(b,0,3,100,17,'Reference-removal experiments\n1,034 development queries × 9 settings\nQueries are reused across settings',fs=8.6)
    b.text(50,48,'Different reference roles; not added to query totals',ha='center',fontsize=7.4)
    c.text(0,100,'(c)  Additional public evaluation',fontweight='bold',fontsize=10,va='top')
    box(c,0,81,100,12,f"Eligible public metadata pool\nn = {n['public_pool']:,}")
    c.text(5,76,'96 used for method confirmation\n1,471 linked assemblies held out',ha='left',va='top',fontsize=8.2)
    arrow(c,(94,81),(94,61))
    box(c,0,49,100,12,f"Public evaluation assemblies\nn = {n['public_evaluation']:,}", '#DCEEF7')
    c.text(5,44,'36 previously exposed during development\nkept separate from source-name comparison',ha='left',va='top',fontsize=7.9)
    arrow(c,(94,49),(94,29))
    box(c,0,17,100,12,f"Source-name comparison\nn = {n['public_source_comparison']:,}", '#DCEEF7')
    c.text(50,12,'100 expert-labelled genomes are a subset of 4,857.\n20 outcome-selected reviews are a subset of 4,821.',ha='center',va='top',fontsize=8)
    fig.text(.5,.118,'Native classifier runs: 946 study + 4,857 public = 5,803 assemblies',ha='center',fontweight='bold',fontsize=9)
    fig.text(.5,.071,'Accession catalogue (Table S01): 1,098 records = 1,074 study + 16 types + 2 outgroups\n'
        '+ 5 initial exclusions + 1 diagnostic addition.\nCatalogue rows are not an analysis denominator.',ha='center',fontsize=8.5,linespacing=1.4)
    fig.text(.5,.032,'Reference-removal settings, expert subsets and fragmented derivatives overlap these datasets.\n'
        'Their classification outcomes are not additional independent genomes.',ha='center',fontsize=8.2,linespacing=1.4)
    save(fig,'figure1',[RUN/'outputs/dataset-flow.json',RUN/'outputs/tables/S01-dataset-flow.tsv'],
         'Study accession selection, phylogenetic subsets, primary/genomovar references and public evaluation. All splits are accession/role accounting; overlapping subsets are explicitly identified.',
         ['Counts reconstructed from saved manifests','No new sequence analysis or label change','Original Figure 1 retained in previous run; category totals remain in main Table 1','Boxes have equal role-based sizes; area does not encode sample size'])
