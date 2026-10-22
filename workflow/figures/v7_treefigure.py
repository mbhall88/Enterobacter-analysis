"""Display saved ANI measurements/memberships and reorient a type-tree figure.

No clustering, ANI measurement or phylogenetic inference is performed.
"""
import collections, csv, itertools, json, math, os, shutil
from pathlib import Path
from prepare import ROOT, RUN, sha
os.environ['MPLCONFIGDIR']=str(RUN/'outputs/mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import Phylo

OUT=RUN/'outputs/figures';DATA=OUT/'source-data'
SOURCES={
 'genomovar_pairs':ROOT/'results/ani-cases/2026-10-02/task12-astra-group-a-populations/inputs/pairs.tsv',
 'genomovar_memberships':ROOT/'results/ani-cases/2026-10-02/task12-astra-genomovars-reporting/outputs/memberships.tsv',
 'subspecies_ranges':ROOT/'results/ani-cases/2026-09-29/task12-astra-kobei-xiangfangensis-review/outputs/subgroup-pairs.tsv',
}
def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_rows(p,data):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
records=[]
def save(fig,name,description,transformations,source_keys):
    for ext in ['png','pdf','svg']:
        options={'metadata':{'CreationDate':None,'ModDate':None}} if ext=='pdf' else ({'metadata':{'Date':None}} if ext=='svg' else {})
        fig.savefig(OUT/f'{name}.{ext}',dpi=300,facecolor='white',**options)
    plt.close(fig)
    records.append({'id':name,'description':description,'transformations':transformations,
                    'source_files':[{'path':str(SOURCES[k]),'sha256':sha(SOURCES[k])} for k in source_keys],
                    'formats':{e:{'path':str(OUT/f'{name}.{e}'),'sha256':sha(OUT/f'{name}.{e}')} for e in ['png','pdf','svg']},
                    'phase':'Internal IJSEM author-review draft; final publication approval pending'})

def tree_figure():
    names=['bac120_fasttree','verticall_fastme','mashtree_consensus']
    display=rows(RUN/'outputs/tables/S02-type-display-map.tsv')
    strains={x['accession']:x['deposited_strain_label'] for x in rows(RUN/'outputs/figures/source-data/figure2-strain-labels.tsv')}
    by={r['accession']:(i+1,r) for i,r in enumerate(display)}
    outgroups={r['accession'] for r in display if 'outgroup' in r['current_project_label']}
    fig=plt.figure(figsize=(10.2,12.5),layout='constrained');grid=fig.add_gridspec(2,3,height_ratios=[1.75,1])
    keyax=fig.add_subplot(grid[1,:]);keyax.axis('off');reports={}
    titles=['(a)  BAC120 / FastTree\nSH-like local support','(b)  Verticall / FastME\nNo branch-support estimates','(c)  Mashtree consensus\nSketch-seed frequency']
    for col,(name,title) in enumerate(zip(names,titles)):
        p=ROOT/f'results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/type-anchor-projection-35/trees/{name}.nwk'
        SOURCES[name]=p;t=Phylo.read(p,'newick');nodes=list(t.find_clades());ix={c:i for i,c in enumerate(nodes)}
        graph={i:{} for i in range(len(nodes))};leaf={ix[c]:c.name for c in t.get_terminals()}
        for c in nodes:
            for child in c.clades:
                graph[ix[c]][ix[child]]=graph[ix[child]][ix[c]]={'length':float(child.branch_length or 0),'support':child.confidence}
        all_tips=set(leaf.values())
        def side(g,start,avoid):
            todo=[start];visited={avoid};found=set()
            while todo:
                n=todo.pop()
                if n in visited:continue
                visited.add(n)
                if n in leaf:found.add(leaf[n])
                todo.extend(v for v in g[n] if v not in visited)
            return found
        def signature(g):
            answer={}
            for a,edges in g.items():
                for b,e in edges.items():
                    if a>=b:continue
                    s=side(g,a,b);other=all_tips-s;key=tuple(sorted(min((s,other),key=lambda v:(len(v),sorted(v)))))
                    record=answer.setdefault(key,{'length':0,'support':set()});record['length']+=e['length']
                    if e['support'] is not None:record['support'].add(float(e['support']))
            return answer
        original=signature(graph)
        m=t.common_ancestor(list(outgroups));assert {c.name for c in m.get_terminals()}==outgroups
        child=ix[m];parent=next(n for n in graph[child] if side(graph,n,child)==all_tips-outgroups)
        edge=graph[parent].pop(child);graph[child].pop(parent);root=len(graph);graph[root]={}
        # Split the existing outgroup stem for display only. Keep its support on
        # one of the two representations so it is not counted twice.
        for node,support in [(child,None),(parent,edge['support'])]:
            e={'length':edge['length']/2,'support':support};graph[root][node]=e;graph[node][root]=e
        final=signature(graph)
        assert original.keys()==final.keys()
        assert all(math.isclose(original[k]['length'],final[k]['length'],rel_tol=1e-12,abs_tol=1e-12) and original[k]['support']==final[k]['support'] for k in original)
        children={};depth={root:0};support={root:None}
        def orient(n,prev):
            children[n]=[x for x in graph[n] if x!=prev]
            children[n].sort(key=lambda x:(side(graph,x,n)==outgroups,sorted(side(graph,x,n))))
            for x in children[n]:
                depth[x]=depth[n]+graph[n][x]['length'];support[x]=graph[n][x]['support'];orient(x,n)
        orient(root,None);ordered=[]
        def tips(n):
            if n in leaf:ordered.append(n)
            for x in children[n]:tips(x)
        tips(root);yy={n:35-i for i,n in enumerate(ordered)}
        def ypos(n):
            if n not in yy:yy[n]=sum(ypos(x) for x in children[n])/len(children[n])
            return yy[n]
        ypos(root)
        # Independent path-length check through the new graph against Bio.Phylo.
        def distance(a,b):
            todo=[(a,None,0)]
            while todo:
                n,prev,d=todo.pop()
                if n==b:return d
                todo.extend((x,n,d+e['length']) for x,e in graph[n].items() if x!=prev)
        by_name={v:k for k,v in leaf.items()};pairs=0
        for a,b in itertools.combinations(sorted(all_tips),2):
            assert math.isclose(distance(by_name[a],by_name[b]),t.distance(a,b),rel_tol=1e-12,abs_tol=1e-12);pairs+=1
        ax=fig.add_subplot(grid[0,col]);maximum=max(depth.values())
        for n in children:
            if children[n]:
                ax.plot([depth[n],depth[n]],[min(yy[x] for x in children[n]),max(yy[x] for x in children[n])],color='#333333',lw=.7)
                for x in children[n]:ax.plot([depth[n],depth[x]],[yy[x],yy[x]],color='#333333',lw=.7)
                if support[n] is not None and support[n]>=.95:ax.plot(depth[n],yy[n],'o',color='#333333',ms=2.4)
            else:
                accession=leaf[n];number,row=by[accession];color='#D55E00' if accession in outgroups else '#333333'
                ax.text(depth[n]+maximum*.025,yy[n],f'T{number:02}'+(' [out]' if accession in outgroups else r'$^{T}$'),fontsize=8.4,color=color,va='center')
        ax.set_ylim(0,36);ax.set_xlim(0,maximum*1.65);ax.set_yticks([]);ax.set_xticks([])
        for s in ax.spines.values():s.set_visible(False)
        ax.set_title(title,loc='left',fontsize=10);length=maximum/5
        ax.plot([0,length],[-.5,-.5],color='#333333',clip_on=False);ax.text(length/2,-1,f'{length:.3f} source distance units',ha='center',fontsize=6.6)
        reports[name]={'input_sha256':sha(p),'tip_count':35,'outgroups':sorted(outgroups),'outgroup_split_present':True,
                       'same_unrooted_splits_lengths_support':True,'pairwise_distance_checks':pairs,
                       'drawing_nodes':[{'id':n,'name':leaf.get(n),'x':depth[n],'y':yy[n],'support':support[n],'children':children[n]} for n in children]}
    for i,r in enumerate(display):
        col=i//18;row=i%18;label=r['current_project_label'].replace('Enterobacter ','E. ')
        aliases={'GCA_001729745.1':' [hoffmannii type]','GCA_040096145.1':' [intestinihominis type]','GCA_052922205.1':' [nototheniae type]'}
        label=label+aliases.get(r['accession'],'')
        keyax.text(col*.5,1-row*.055,f'T{i+1:02}'+(' [out]' if r['accession'] in outgroups else r'$^{T}$')+'  '+label+'\n       '+strains[r['accession']]+'; '+r['accession'],fontsize=6.4,va='top',transform=keyax.transAxes)
    fig.suptitle('Type-strain phylogenies oriented on the shared outgroup edge\nDots: source local support or sketch frequency ≥0.95',fontsize=10)
    (OUT/'figure2-oriented-geometry.json').write_text(json.dumps(reports,indent=2)+'\n')
    save(fig,'figure2','Existing 35-tip projections displayed with the same outgroup orientation',
         ['Place display root halfway along the existing edge separating the two outgroups',
          'Retain every unrooted split, path length and edge support; original Newick files unchanged',
          'This 35-tip display does not infer a root for all 1,089 tips'],names)
    return reports


if __name__=='__main__':
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'Enterobacter-v7-types'}):
        trees=tree_figure()
    (RUN/'checks/v7-figure2.json').write_text(json.dumps({'trees':trees,'distance_checks':sum(x['pairwise_distance_checks'] for x in trees.values()),'all_geometry_preserved':all(x['same_unrooted_splits_lengths_support'] for x in trees.values()),'label_change':'Superscript T marks selected historical type genomes; key includes every assembly accession. Tree data and scientific labels unchanged.'},indent=2)+'\n')
    print('All 1,785 pairwise path lengths, unrooted splits and support retained; 35 accessions in figure key.')
