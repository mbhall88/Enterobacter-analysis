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
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    records.append({'id':name,'description':description,'transformations':transformations,
                    'source_files':[{'path':str(SOURCES[k]),'sha256':sha(SOURCES[k])} for k in source_keys],
                    'formats':{e:{'path':str(OUT/f'{name}.{e}'),'sha256':sha(OUT/f'{name}.{e}')} for e in ['png','pdf','svg']},
                    'phase':'Internal IJSEM author-review draft; final publication approval pending'})

def genomovars():
    members=sorted(rows(SOURCES['genomovar_memberships']),key=lambda r:(r['genomovar_id'],r['accession']))
    pairs=rows(SOURCES['genomovar_pairs']);assert len(members)==173 and len(pairs)==14878
    ids=[r['accession'] for r in members];ix={a:i for i,a in enumerate(ids)}
    groups=[r['genomovar_id'] for r in members];sizes=collections.Counter(groups)
    matrix=np.full((173,173),np.nan);cross=[];within=[]
    for r in pairs:
        i,j=ix[r['a']],ix[r['b']];assert np.isnan(matrix[i,j]) and i!=j
        v=float(r['ani']);matrix[i,j]=matrix[j,i]=v
        assert min(float(r['af_a']),float(r['af_b']))>=80
        if groups[i]==groups[j]:within.append(v)
        elif v>=99.5:cross.append((i,j,v))
    assert len(sizes)==25 and min(within)>=99.5 and len(cross)==11
    assert np.isnan(matrix).sum()==173
    write_rows(DATA/'figureS9-genome-order.tsv',[{'matrix_index':i+1,**r} for i,r in enumerate(members)])
    shutil.copy2(SOURCES['genomovar_pairs'],DATA/'figureS9-pairs.tsv')
    write_rows(DATA/'figureS9-above-threshold-cross-group.tsv',[
        {'a':ids[i],'b':ids[j],'genomovar_a':groups[i],'genomovar_b':groups[j],'ani':v} for i,j,v in cross])
    fig,(a,b)=plt.subplots(1,2,figsize=(9.2,5.7),layout='constrained',gridspec_kw={'width_ratios':[1.55,1]})
    cmap=matplotlib.colormaps['viridis'].with_extremes(bad='#eeeeee')
    im=a.imshow(matrix,vmin=98,vmax=100,cmap=cmap,interpolation='nearest',origin='upper')
    assert np.nanmin(matrix)>=98
    offset=0
    for group,size in sizes.items():
        a.add_patch(Rectangle((offset-.5,offset-.5),size,size,fill=False,edgecolor='#ffffff',lw=.8))
        offset+=size
    for i,j,v in cross:
        a.plot(j,i,marker='o',mfc='none',mec='#D55E00',ms=4,lw=0,mew=.8)
        a.plot(i,j,marker='o',mfc='none',mec='#D55E00',ms=4,lw=0,mew=.8)
    ticks=[0,49,99,149,172];a.set_xticks(ticks,[1,50,100,150,173]);a.set_yticks(ticks,[1,50,100,150,173])
    a.set_xlabel('Genome index, ordered by genomovar');a.set_ylabel('Genome index, ordered by genomovar')
    a.set_title('a  ANI among the 173 definition genomes',loc='left',fontsize=10)
    fig.colorbar(im,ax=a,orientation='horizontal',fraction=.045,pad=.14,label='skani ANI (%)',ticks=[98,98.5,99,99.5,100])
    labels=list(sizes);counts=list(sizes.values());ys=np.arange(25)
    b.barh(ys,counts,color='#0072B2',edgecolor='#333333',height=.7)
    b.set_yticks(ys,labels,fontsize=7);b.invert_yaxis();b.set_xlim(0,104)
    b.set_xlabel('Genomes');b.set_title('b  Membership of the 25 genomovars',loc='left',fontsize=10)
    for y,n in zip(ys,counts):b.text(n+1,y,str(n),va='center',fontsize=7)
    fig.suptitle('Genomovars within the proposed E. xiangfangensis subsp. hoffmannii',fontsize=11)
    save(fig,'figureS9','ANI heatmap and counts for the unchanged 25-genomovar partition',
         ['Genomes ordered by frozen genomovar ID then accession; no reclustering',
          'All 14,878 measured pairs shown; unmeasured diagonal is grey',
          'White boxes show within-group blocks; orange circles mark the eleven between-group pairs at or above 99.5%',
          'Fixed colour limits 98–100%; all measured values lie within this range'],
         ['genomovar_pairs','genomovar_memberships'])
    return {'genomes':len(ids),'pairs':len(pairs),'groups':len(sizes),'singletons':sum(n==1 for n in counts),'largest':max(counts),'within_min':min(within),'between_at_or_above':len(cross),'minimum_ani':float(np.nanmin(matrix))}

def subspecies():
    data=[r for r in rows(SOURCES['subspecies_ranges']) if r['scope']=='without_peripheral437']
    assert len(data)==6
    data.sort(key=lambda r:(r['group_a']!=r['group_b'],r['group_a'],r['group_b']))
    write_rows(DATA/'figureS10-subspecies-ranges.tsv',data)
    fig,ax=plt.subplots(figsize=(8.5,4.7),layout='constrained')
    for i,r in enumerate(data):
        same=r['group_a']==r['group_b'];color='#0072B2' if same else '#333333'
        low,median,high=[float(r[x]) for x in ['min_ani','median_ani','max_ani']]
        ax.plot([low,high],[i,i],color=color,lw=2)
        ax.plot([low,high],[i,i],linestyle='none',marker='|',color=color,ms=10)
        ax.plot(median,i,marker='o' if same else 'D',color=color,ms=5)
        ax.text(100.08,i,f"n = {int(r['pairs']):,}",fontsize=8,va='center')
    labels=[r['group_a'] if r['group_a']==r['group_b'] else r['group_a']+' / '+r['group_b'] for r in data]
    ax.set_yticks(range(6),labels);ax.invert_yaxis();ax.set_xlim(96.5,100.8)
    ax.set_xticks([96.5,97,97.5,98,98.5,99,99.5,100]);ax.set_xlabel('skani ANI (%)')
    ax.axhline(2.5,color='#999999',lw=.6,ls='--')
    ax.set_title('Within-subspecies (circles) and between-subspecies (diamonds) ANI\nRanges and medians among 437 assigned genomes',fontsize=11)
    save(fig,'figureS10','Measured minimum, median and maximum ANI within and between three assigned subspecies',
         ['Use only the saved without_peripheral437 rows: 18 oharae, 256 steigerwaltii, 163 xiangfangensis',
          'The species-only peripheral genome and its nearest-type sensitivity bin are excluded',
          'Lines show observed min–max and points show medians, not confidence intervals',
          'Pairwise observations share genomes; no statistical independence assumed'],['subspecies_ranges'])
    return {'rows':len(data),'pairs':sum(int(r['pairs']) for r in data),'genomes':437}

def tree_figure():
    names=['bac120_fasttree','verticall_fastme','mashtree_consensus']
    display=rows(RUN/'outputs/tables/S02-type-display-map.tsv')
    by={r['accession']:(i+1,r) for i,r in enumerate(display)}
    outgroups={r['accession'] for r in display if 'outgroup' in r['current_project_label']}
    fig=plt.figure(figsize=(8,10.4),layout='constrained');grid=fig.add_gridspec(2,3,height_ratios=[2.3,1])
    keyax=fig.add_subplot(grid[1,:]);keyax.axis('off');reports={}
    titles=['a  BAC120 / FastTree\nSH-like local support','b  Verticall / FastME\nNo branch-support estimates','c  Mashtree consensus\nSketch-seed frequency']
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
                ax.text(depth[n]+maximum*.025,yy[n],f'T{number:02}'+(' [out]' if accession in outgroups else ''),fontsize=8.4,color=color,va='center')
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
        label=label.replace('E. xiangfangensis subsp. hoffmannii','E. x. subsp. hoffmannii')+aliases.get(r['accession'],'')
        keyax.text(col*.5,1-row*.055,f'T{i+1:02}  '+label,fontsize=6.7,va='top',transform=keyax.transAxes)
    fig.suptitle('Type-strain phylogenies oriented on the shared outgroup edge\nDots: source local support or sketch frequency ≥0.95',fontsize=10)
    (OUT/'figure2-oriented-geometry.json').write_text(json.dumps(reports,indent=2)+'\n')
    save(fig,'figure2','Existing 35-tip projections displayed with the same outgroup orientation',
         ['Place display root halfway along the existing edge separating the two outgroups',
          'Retain every unrooted split, path length and edge support; original Newick files unchanged',
          'This 35-tip display does not infer a root for all 1,089 tips'],names)
    return reports

def main():
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'}):
        gv=genomovars();sub=subspecies();trees=tree_figure()
    (OUT/'v5-figure-manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    (RUN/'checks/v5-figure-data.json').write_text(json.dumps({'genomovars':gv,'subspecies':sub,
        'tree_distance_checks':sum(r['pairwise_distance_checks'] for r in trees.values()),
        'all_tree_splits_lengths_support_preserved':all(r['same_unrooted_splits_lengths_support'] for r in trees.values()),
        'matplotlib':matplotlib.__version__},indent=2)+'\n')
    (RUN/'inputs/v5-figure-sources.json').write_text(json.dumps({k:{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for k,p in SOURCES.items()},indent=2)+'\n')
    print(json.dumps({'genomovars':gv,'subspecies':sub,'tree_pairwise_checks':1785},indent=2))
if __name__=='__main__':main()
