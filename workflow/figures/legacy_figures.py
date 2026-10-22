"""Render saved results. No inference, alignment, sequence input or source edits."""
import csv, json, os
from pathlib import Path
from prepare import RUN, ROOT
os.environ['MPLCONFIGDIR']=str(RUN/'outputs/mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib.patches import FancyBboxPatch
from prepare import sha
OUT=RUN/'outputs/figures';OUT.mkdir(exist_ok=True)
TAB=RUN/'outputs/tables'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white'})
records=[]
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def save(fig,name,sources,description,transforms):
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    records.append({'id':name,'sources':[{'path':str(p),'sha256':sha(p)} for p in sources],'description':description,'transformations':transforms,'formats':['pdf','svg','png'],'uncertainty':'No new interval or inference. Retained counts/point measurements only.','publisher_phase':'IJSEM author-review draft; final author/editorial figure approval pending'})
def figure1():
    from flow_figure import draw
    draw(save,RUN)
def figure2():
    mp=rows(TAB/'S02-type-display-map.tsv');by={r['accession']:(i+1,r) for i,r in enumerate(mp)}
    fig=plt.figure(figsize=(8.0,10.4),layout='constrained')
    grid=fig.add_gridspec(2,3,height_ratios=[2.3,1.0])
    axes=[fig.add_subplot(grid[0,i]) for i in range(3)]
    keyax=fig.add_subplot(grid[1,:]);keyax.axis('off')
    methods=['bac120_fasttree','verticall_fastme','mashtree_consensus']
    titles=['A  BAC120 / FastTree\nSH-like local support','B  Verticall / FastME\nNo branch-support estimates','C  Mashtree consensus\nSketch-seed frequency']
    source=[];geometry={}
    for ax,name,title in zip(axes,methods,titles):
        p=ROOT/f'results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/type-anchor-projection-35/trees/{name}.nwk'
        tree=Phylo.read(p,'newick');terms=tree.get_terminals();assert {c.name for c in terms}==set(by)
        dep=tree.depths();yy={c:len(terms)-i for i,c in enumerate(terms)}
        def assign(c):
            if c not in yy:yy[c]=sum(assign(ch) for ch in c.clades)/len(c.clades)
            return yy[c]
        assign(tree.root)
        maximum=max(dep.values());geometry[name]={'input_sha256':sha(p),'tips':sorted(by),'branches':[{'name':c.name,'length':c.branch_length,'confidence':c.confidence,'x':dep[c],'y':yy[c]} for c in tree.find_clades()]}
        def draw(c):
            if c.clades:
                ax.plot([dep[c],dep[c]],[min(yy[ch] for ch in c.clades),max(yy[ch] for ch in c.clades)],color='#333333',lw=.7)
                for ch in c.clades:
                    ax.plot([dep[c],dep[ch]],[yy[ch],yy[ch]],color='#333333',lw=.7);draw(ch)
                if c.confidence is not None and c.confidence>=.95:
                    ax.plot(dep[c],yy[c],marker='o',color='#333333',ms=2.4)
            else:
                no,r=by[c.name];label=r['current_project_label'].replace('Enterobacter ','').replace(' subsp. ',' s. ')
                label=label.replace('xiangfangensis s. ','xian. s. ').replace('roggenkampii s. ','rogg. s. ').replace('chengduensis s. ','cheng. s. ')
                aliases={'GCA_001729745.1':'hoffmannii type','GCA_040096145.1':'intestinihominis type','GCA_052922205.1':'nototheniae type'}
                if c.name in aliases: label=aliases[c.name]+' → xian. s. hoffmannii'
                label=label.replace('Cronobacter sakazakii (outgroup)','C. sakazakii [out]').replace('Klebsiella aerogenes (outgroup)','K. aerogenes [out]')
                color='#0072B2' if c.name in ['GCA_001729745.1','GCA_040096145.1','GCA_052922205.1'] else '#333333'
                ax.text(dep[c]+maximum*.025,yy[c],f'T{no:02}',fontsize=9,va='center',color=color)
        draw(tree.root);ax.set_ylim(0,36);ax.set_xlim(0,maximum*1.3);ax.set_yticks([]);ax.set_title(title,fontsize=10,loc='left')
        ax.spines['left'].set_visible(False);ax.spines['bottom'].set_visible(False);ax.set_xticks([])
        length=maximum/5;ax.plot([0,length],[-.5,-.5],color='#333333',clip_on=False)
        ax.text(length/2,-1.0,f'{length:.3f} source distance units',ha='center',fontsize=6.6)
        source.append(p)
    for i,r in enumerate(mp):
        col=i//18;row=i%18
        label=r['current_project_label'].replace('Enterobacter ','E. ')
        aliases={'GCA_001729745.1':' [hoffmannii type]', 'GCA_040096145.1':' [intestinihominis type]', 'GCA_052922205.1':' [nototheniae type]'}
        label=label.replace('E. xiangfangensis subsp. hoffmannii','E. x. subsp. hoffmannii')+aliases.get(r['accession'],'')
        keyax.text(col*.50,1-row*.055,f'T{i+1:02}  '+label,fontsize=6.7,va='top',transform=keyax.transAxes)
    fig.suptitle('Type-strain phylogenies with revised classifications\nDots: source local support or sketch frequency ≥0.95; original orientation retained',fontsize=10)
    (OUT/'figure2-geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
    write=[{'type_display_id':f'T{i+1:02}',**r} for i,r in enumerate(mp)]
    with (TAB/'S02-type-display-map.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(write[0]),delimiter='\t');w.writeheader();w.writerows(write)
    save(fig,'figure2',source+[TAB/'S02-type-display-map.tsv'],'Three unchanged 35-tip type projections. Abbreviated current names and unique T IDs link to full historical/current map. Dots have different method semantics; proposed names remain draft.',['Tip IDs and full name key added; panel layout revised for legibility','Input branch lengths, tip sets and confidence fields retained','Arbitrary original orientation; no new rooting claim'])
def figure3():
    ds=json.loads((TAB/'S08-dataset-setting-summary.json').read_text())
    arm=[r for r in ds if r['dataset']=='layer2' and r['applicable'] and r['setting'].endswith('omit_reviewed_subspecies')]
    fig,(a,b)=plt.subplots(1,2,figsize=(10.3,5.2),layout='constrained',gridspec_kw={'width_ratios':[.9,1.5]})
    labels=[{'self_alias':'Query removed','clone_block':'Clone group removed','split_block':'Broader group removed'}.get(r['setting'].split('|')[0],r['setting'].split('|')[0].replace('_',' ')) for r in arm]
    rescue=[r['counts']['rescued'] for r in arm];veto=[r['counts']['contexts_with_method_veto'] for r in arm]
    a.barh(labels,rescue,color='#0072B2',label='Species assigned with ANIm',edgecolor='#333333')
    a.barh(labels,veto,left=rescue,color='#D55E00',hatch='///',label='Unassigned: methods disagree',edgecolor='#333333')
    for i,(n,v) in enumerate(zip(rescue,veto)):
        a.text(n/2,i,str(n),color='white',ha='center',va='center')
        if v:a.text(n+v/2 if v>=8 else n+v+2,i,str(v),color='black',ha='center' if v>=8 else 'left',va='center')
    a.set_xlabel('Classifications (queries reused across settings)');a.set_xlim(0,90);a.set_title('A  Subspecies references removed',loc='left',fontsize=11,fontweight='bold')
    a.legend(loc='lower left',bbox_to_anchor=(-.02,-.29),fontsize=8)
    cat=json.loads((TAB/'S08-genomovar-category-counts.json').read_text())
    mappings=[
      ('Unique group match','gv_source_group','gv_assigned','#0072B2',''),
      ('Singleton match',None,'gv_assigned_singleton','#009E73','++'),
      ('Unassessable singleton','gv_unassessable_singleton',None,'#FFFFFF','..'),
      ('Multiple groups match','gv_withheld_overlap','gv_withheld_overlap','#D55E00','///'),
      ('Missing comparisons','gv_withheld_incomplete',None,'#BBBBBB','xx'),
      ('Insufficient member coverage',None,'gv_withheld_member_af','#CC79A7','--'),
      ('No compatible group',None,'gv_no_compatible','#777777','\\\\'),
    ]
    datasets=['GV1b','GV2','GV3-public'];ys=list(range(3));left=[0.,0.,0.]
    for label,key,pubkey,col,hatch in mappings:
        nums=[cat[d].get(pubkey if d=='GV3-public' else key,0) if (pubkey if d=='GV3-public' else key) else 0 for d in datasets]
        totals=[173,173,749];perc=[100*n/t for n,t in zip(nums,totals)]
        b.barh(ys,perc,left=left,color=col,edgecolor='#333333',hatch=hatch,label=label)
        left=[x+y for x,y in zip(left,perc)]
    assert all(abs(x-100)<1e-8 for x in left)
    b.set_yticks(ys,['Leave-one-out\nn=173','Related genomes removed\nn=173','Public assemblies\nn=749']);b.invert_yaxis();b.set_xlim(0,100);b.set_xlabel('Classification outcomes (%)');b.set_title('B  Genomovar matching',loc='left',fontsize=11,fontweight='bold')
    b.legend(loc='upper left',bbox_to_anchor=(0,-.18),ncol=2,fontsize=8)
    save(fig,'figure3',[TAB/'S08-dataset-setting-summary.json',TAB/'S08-genomovar-category-counts.json'],'Only three overlapping development arms have rescues. Separate genomovar strata show unique, singleton, overlap, missing, coverage and no-compatible outcomes.',['Counts from accepted reporting correction','Bars normalized separately by 173, 173 and 749 records','No pooled inference or independent confidence interval'])
def figure4():
    data=[r for r in rows(TAB/'S09-Kraken-metrics.tsv') if r['dimension']=='resolution']
    cats=['species_or_below','below_genus_without_named_species','genus_only','above_genus_or_unresolved','unclassified']
    colors=['#0072B2','#56B4E9','#BBBBBB','#E69F00','#FFFFFF'];hatches=['','..','///','xx','++']
    settings=[('native','0'),('native','0.1'),('fragments150','0'),('fragments150','0.1')]
    fig,ax=plt.subplots(figsize=(8.4,4.4),layout='constrained');left=[0.]*4
    for cat,col,hatch in zip(cats,colors,hatches):
        # Sparse metrics omit zero rows. Check the full category sum below.
        vals=[sum(float(r['percent_represented_bases']) for r in data if r['input_mode']==m and r['confidence']==c and r['category']==cat) for m,c in settings]
        ax.barh(range(4),vals,left=left,label=cat.replace('_',' '),color=col,hatch=hatch,edgecolor='#333333')
        if cat=='species_or_below':
            for i,v in enumerate(vals):ax.text(v/2,i,f'{v:.3f}%',color='white',ha='center',va='center',fontsize=9)
        left=[a+b for a,b in zip(left,vals)]
    assert all(abs(x-100)<1e-6 for x in left)
    ax.set_yticks(range(4),['Native; confidence 0','Native; confidence 0.1','≤150-bp fragments; confidence 0','≤150-bp fragments; confidence 0.1']);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Input bases by assigned rank (%)')
    ax.legend(loc='upper left',bbox_to_anchor=(0,-.18),ncol=2,fontsize=8);ax.set_title('Kraken2 taxonomic resolution by sequence length and confidence',loc='left',fontweight='bold',fontsize=11)
    save(fig,'figure4',[TAB/'S09-Kraken-metrics.tsv'],'Four rank-retention profiles over identical base totals. Artificial fragments and native contigs differ in species/genus representation. This is not isolate accuracy.',['Copied accepted exact metrics','Representation/confidence displayed as separate settings','Full sequence-length weighting; no basewise-support inference'])
if __name__=='__main__':
    figure1();figure2();figure3();figure4()
    (OUT/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    (OUT/'environment.json').write_text(json.dumps({'matplotlib':matplotlib.__version__,'backend':matplotlib.get_backend()},indent=2)+'\n')

    # Figure2 adds stable display IDs; bind the final table bytes after that step.
    (RUN/'outputs/table-manifest.json').write_text(json.dumps([{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(TAB.iterdir()) if p.is_file()],indent=2)+'\n')
