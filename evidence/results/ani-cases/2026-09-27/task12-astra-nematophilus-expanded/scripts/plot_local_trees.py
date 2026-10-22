"""Draw two reference sensitivities with every diagnostic tip visible."""
from pathlib import Path
import argparse, csv, json, itertools, math
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p=argparse.ArgumentParser();p.add_argument('--condition',choices=['primary','masked','RY'],required=True);a=p.parse_args()
r=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((r/'inputs/panel-qc.tsv').open(),delimiter='\t'))
rows=[x for x in rows if x['local_primary']=='True'];ids={x['accession'] for x in rows}
group={x['accession']:('nematophilus' if x['accession']=='GCA_026344075.1' or x['group']=='nematophilus_candidate' else x['group']) for x in rows}
groups={g:{q for q in ids if group[q]==g} for g in set(group.values())}
colours={'bugandensis':'#0072B2','nematophilus':'#D55E00','boundary':'#8B51A2','context':'#707070'}
names={'GCA_019046905.1':'bugandensis TYPE','GCA_026344075.1':'nematophilus TYPE','GCA_043966805.1':'ECC-219','GCA_030510045.1':'strain 93 [new]','GCA_001562175.1':'MBRL1077 [new]','GCA_024160165.1':'9i [new]','GCA_034107265.1':'168J2 [new]','GCA_047525435.1':'4EL [new]','GCA_055683465.1':'WI-CRE-C84 [new]','GCA_003594915.1':'chuandaensis TYPE'}
def canonical(s):
    other=ids-s
    return tuple(sorted(s if len(s)<len(other) else other if len(s)>len(other) else min(tuple(sorted(s)),tuple(sorted(other)))))
fig,axes=plt.subplots(1,2,figsize=(20,13));checks=[]
for ax,ref in zip(axes,['GCA_019046905.1','GCA_026344075.1']):
    o=r/'local'/ref
    path=o/('tree/primary.treefile' if a.condition=='primary' else 'tests/'+a.condition+'.treefile')
    tree=Phylo.read(path,'newick');assert {x.name for x in tree.get_terminals()}==ids
    # Capture labels by unrooted split before changing the display root.
    support={canonical({x.name for x in c.get_terminals()}):(c.name or str(c.confidence or '')) for c in tree.get_nonterminals()}
    before={(x,y):tree.distance(x,y) for x,y in itertools.combinations(sorted(ids),2)}
    context=tree.common_ancestor(sorted(groups['context']))
    if {x.name for x in context.get_terminals()}!=groups['context']:
        # Root first on a focal tip so the context clade is represented as descendants.
        tree.root_with_outgroup('GCA_019046905.1')
        context=tree.common_ancestor(sorted(groups['context']))
    assert {x.name for x in context.get_terminals()}==groups['context']
    tree.root_with_outgroup(context);tree.ladderize(reverse=True)
    assert all(math.isclose(tree.distance(x,y),d,rel_tol=1e-10,abs_tol=1e-12) for (x,y),d in before.items())
    tips=tree.get_terminals();ys={c:i for i,c in enumerate(tips)};xs=tree.depths()
    for c in tree.get_nonterminals(order='postorder'):ys[c]=sum(ys[k] for k in c.clades)/len(c.clades)
    maximum=max(xs.values());pad=maximum*.025
    wanted=[groups[g] for g in groups]+[groups['boundary']|groups['nematophilus'],groups['boundary']|groups['bugandensis'],groups['nematophilus']|groups['bugandensis']]
    labelled=set()
    for c in tree.find_clades():
        members={x.name for x in c.get_terminals()};gs={group[q] for q in members};colour=colours[next(iter(gs))] if len(gs)==1 else '#222222'
        for child in c.clades:ax.plot([xs[c],xs[child]],[ys[child]]*2,color=colour,lw=1.2)
        if c.clades:ax.plot([xs[c]]*2,[min(ys[k] for k in c.clades),max(ys[k] for k in c.clades)],color=colour,lw=1.2)
        key=canonical(members)
        if c.clades and key in [canonical(v) for v in wanted] and key not in labelled:
            label=support.get(key,'');labelled.add(key)
            if label:ax.text(xs[c],ys[c]-.35,label,fontsize=8,ha='right',va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=.5))
        if not c.clades:
            ax.scatter([xs[c]],[ys[c]],s=12,color=colour,zorder=4)
            text=c.name+('  '+names[c.name] if c.name in names else '')
            ax.text(xs[c]+pad,ys[c],text,fontsize=8,color=colour,va='center',weight='bold' if c.name in names else 'normal')
    validation=json.loads((o/'validation.json').read_text())
    bp=validation['strict_core_bases'] if a.condition!='masked' else json.loads((o/'tests/recombination-summary.json').read_text())['retained_core_bases']
    ax.set_title(names[ref]+' reference\n'+f'{bp:,} shared sites; 38 genomes',fontsize=13)
    ax.set_ylim(38,-1);ax.set_xlim(-pad,maximum*2.3);ax.set_yticks([])
    ax.set_xlabel('Substitutions / site'+(' (binary R/Y model)' if a.condition=='RY' else ''))
    for side in ['left','right','top']:ax.spines[side].set_visible(False)
    checks.append(dict(reference=ref,tree=str(path),tips=38,patristic_checks=len(before),condition=a.condition))
title={'primary':'Local nucleotide phylogenies: unmasked shared core','masked':'Local nucleotide phylogenies: recombination-masked shared core','RY':'Composition sensitivity: R/Y-recoded shared core'}[a.condition]
fig.suptitle(title,fontsize=18,y=.99)
fig.legend(handles=[Line2D([0],[0],color=c,lw=3,label=f'{g} ({len(groups[g])})') for g,c in colours.items()],loc='lower center',ncol=4,bbox_to_anchor=(.5,.038))
fig.text(.5,.018,'Selected branches: SH-aLRT / ultrafast bootstrap (%). Rooted on the four neighbouring genomes for display.\nColour denotes diagnostic group membership, not a final species assignment; no tips collapsed.',ha='center',fontsize=10)
fig.subplots_adjust(top=.92,bottom=.10,left=.035,right=.99,wspace=.07)
for ext in ['png','pdf']:fig.savefig(r/f'figures/local-trees-{a.condition}.{ext}',dpi=180)
(r/f'figures/local-trees-{a.condition}-validation.json').write_text(json.dumps(dict(status='PASS',checks=checks),indent=2)+'\n')
