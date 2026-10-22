"""Author-facing display of frozen comparisons; no inference or reclassification."""
from collections import Counter
import json, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from prepare import RUN, ROOT, rows, sha
from v6_data import kraken_data, write_tsv, SCORING

OUT=RUN/'outputs/figures'; DATA=OUT/'source-data'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                    'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'Enterobacter-v6'})

def save(fig,stem):
    for ext in ['png','pdf','svg']:
        kwargs={'dpi':300,'bbox_inches':'tight','facecolor':'white'}
        if ext=='pdf':kwargs['metadata']={'CreationDate':None,'ModDate':None}
        if ext=='svg':kwargs['metadata']={'Date':None}
        fig.savefig(OUT/(stem+'.'+ext),**kwargs)
    plt.close(fig)

def short(s):return s.replace('Enterobacter ','E. ').replace('Lelliottia ','L. ')
def label(s):return short(s).replace(' subsp. ','\nsubsp. ')

def matrix(data,projects,names):
    cells=Counter({(r['study_assignment'],r['kraken_name']):r['bases'] for r in data})
    totals=Counter()
    for r in data:totals[r['study_assignment']]+=r['bases']
    counts=np.zeros((len(projects),len(names)+1),dtype=np.int64)
    for i,p in enumerate(projects):
        for j,n in enumerate(names):counts[i,j]=cells[p,n]
        counts[i,-1]=totals[p]-sum(counts[i,:-1])
    assert (counts>=0).all() and list(counts.sum(axis=1))==[totals[p] for p in projects]
    percent=counts/np.array([totals[p] for p in projects])[:,None]*100
    return counts,percent,totals

def heatmap(ax,data,projects,names,full=False):
    counts,percent,totals=matrix(data,projects,names)
    im=ax.imshow(percent,cmap='Blues',vmin=0,vmax=100,aspect='auto',interpolation='none')
    ax.set_xticks(range(len(names)+1),[label(n) for n in names]+['All other\ncalls'],rotation=55,ha='right',fontsize=8 if full else 9)
    ax.set_yticks(range(len(projects)),[label(p) for p in projects],fontsize=8 if full else 10)
    ax.set_xticks(np.arange(-.5,len(names)+1),minor=True);ax.set_yticks(np.arange(-.5,len(projects)),minor=True)
    ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',bottom=False,left=False)
    for i in range(len(projects)):
        for j in range(len(names)+1):
            value=percent[i,j]
            if value>=1:
                ax.text(j,i,f'{value:.0f}' if full else f'{value:.1f}',ha='center',va='center',fontsize=6 if full else 9,color='white' if value>55 else '#222222')
    ax.set_xlabel('Name returned by Kraken2',labelpad=10)
    ax.set_ylabel('Assignment in this study',labelpad=12)
    return im,counts,percent

def kraken():
    data,transformation=kraken_data()
    write_tsv(DATA/'figure3-and-S13-kraken-taxonomy.tsv',data)
    write_tsv(DATA/'figure3-and-S13-taxonomy-transformation.tsv',transformation)
    projects=['Enterobacter hormaechei']+['Enterobacter xiangfangensis subsp. '+s for s in ['hoffmannii','oharae','steigerwaltii','xiangfangensis']]+['Enterobacter xiangfangensis']
    names=['Enterobacter hormaechei subsp. '+s for s in ['hormaechei','hoffmannii','oharae','steigerwaltii','xiangfangensis']]+['Enterobacter hormaechei','Enterobacter xiangfangensis']
    fig,ax=plt.subplots(figsize=(8.5,5.7));fig.subplots_adjust(left=.31,right=.91,bottom=.35,top=.9)
    im,counts,percent=heatmap(ax,data,projects,names)
    ax.set_title('Kraken names in the hormaechei-related lineages',fontsize=12,pad=12)
    cb=fig.colorbar(im,ax=ax,location='right',fraction=.04,pad=.035);cb.set_label('Input bases within study group (%)',fontsize=9)
    save(fig,'figure3')
    display=[]
    for i,p in enumerate(projects):
        for j,n in enumerate(names+['All other calls']):display.append({'study_assignment':p,'kraken_name':n,'bases':int(counts[i,j]),'percent_within_study_group':float(percent[i,j])})
    write_tsv(DATA/'figure3-displayed-cells.tsv',display)
    # Show every study category. Retain its dominant call and every named column
    # from the main panel; the residual column preserves all lesser calls.
    allprojects=sorted({r['study_assignment'] for r in data})
    dominant={max([r for r in data if r['study_assignment']==p],key=lambda r:(r['bases'],r['kraken_name']))['kraken_name'] for p in allprojects}
    allnames=sorted(dominant | set(names))
    fig,ax=plt.subplots(figsize=(14.5,13.5));fig.subplots_adjust(left=.25,right=.95,bottom=.23,top=.95)
    im,counts,percent=heatmap(ax,data,allprojects,allnames,True)
    ax.set_title('Kraken names versus study assignments: all categories',fontsize=14,pad=12)
    cb=fig.colorbar(im,ax=ax,fraction=.018,pad=.025);cb.set_label('Input bases within study group (%)',fontsize=10)
    save(fig,'figureS13')
    full=[]
    for i,p in enumerate(allprojects):
        for j,n in enumerate(allnames+['All other calls']):full.append({'study_assignment':p,'kraken_name':n,'bases':int(counts[i,j]),'percent_within_study_group':float(percent[i,j])})
    write_tsv(DATA/'figureS13-displayed-cells.tsv',full)
    return {'native_bases':sum(r['bases'] for r in data),'all_study_categories':len(allprojects),'all_raw_ancestor_columns':len({r['kraken_name'] for r in data}),'displayed_named_columns':len(allnames),'main_study_categories':projects,'column_selection':'S13: dominant call in any study row plus all named Figure 3 columns; all remainder retained in residual column; Figure 3: explicit hormaechei-related names; other calls retained'}

def pterobacter():
    study=rows(SCORING/'layer1/scored-genomes.tsv');public=rows(SCORING/'layer3/scored-genomes.tsv')
    named=[r for r in study if r['species_outcome'] in ('correct','withheld')]
    reviewed=[r for r in named if r['label_status']=='approved_review_treatment']
    retained=[r for r in named if r['label_status']=='source_name_retained_not_individually_adjudicated']
    scored=[r for r in public if r['species_outcome'] in ('consistent','discordant','withheld')]
    groups=[('Reviewed study assignments',reviewed),('Retained study names',retained),('Public names',scored)]
    agg=[]
    for name,data in groups:
        c=Counter(r['species_outcome'] for r in data)
        agg.append({'comparison':name,'rank':'species','denominator':len(data),'matching':c['correct']+c['consistent'],'different':c['discordant'],'species_only':0,'withheld':c['withheld']})
    sub=[r for r in study if r['subspecies_outcome'] in ('correct','species_only','withheld_at_species')];c=Counter(r['subspecies_outcome'] for r in sub)
    agg.append({'comparison':'Study subspecies assignments','rank':'subspecies','denominator':len(sub),'matching':c['correct'],'different':0,'species_only':c['species_only'],'withheld':c['withheld_at_species']})
    assert [(r['denominator'],r['matching'],r['different'],r['species_only'],r['withheld']) for r in agg]==[(836,828,0,0,8),(86,86,0,0,0),(4387,4277,56,0,54),(689,677,0,9,3)]
    write_tsv(DATA/'figure4-classification-agreement.tsv',agg)
    colors={'matching':'#009E73','different':'#D55E00','species_only':'#E8BC66','withheld':'#92979D'}
    fig,(ax,bx)=plt.subplots(2,1,figsize=(8.8,5.1),gridspec_kw={'height_ratios':[3,1.7]})
    fig.subplots_adjust(left=.28,right=.78,top=.87,bottom=.13,hspace=.9)
    for a,records,title in [(ax,agg[:3],'a  Species agreement'),(bx,agg[3:],'b  Subspecies agreement')]:
        for i,r in enumerate(records):
            left=0
            for cat in colors:
                value=r[cat]/r['denominator']*100
                a.barh(i,value,left=left,color=colors[cat],height=.55,label=cat if i==0 else None)
                left+=value
            a.text(50,i,f"{r['matching']:,} matching",color='white',ha='center',va='center',fontsize=10)
            other=[]
            if r['different']:other.append(f"{r['different']} different")
            if r['species_only']:other.append(f"{r['species_only']} species only")
            other.append(f"{r['withheld']} withheld")
            a.text(102,i,'\n'.join(other),ha='left',va='center',fontsize=9)
        a.set_yticks(range(len(records)),[f"{r['comparison']}\n(n = {r['denominator']:,})" for r in records])
        a.set_xlim(0,100);a.set_xticks([0,25,50,75,100]);a.invert_yaxis()
        a.set_xlabel('Genomes with a comparable label (%)');a.set_title(title,loc='left',fontweight='bold',fontsize=12,pad=12)
        a.spines['left'].set_visible(False);a.tick_params(axis='y',length=0)
    handles=[plt.Rectangle((0,0),1,1,color=colors[k]) for k in colors]
    fig.legend(handles,['Matching','Different','Species only','Withheld'],ncol=4,loc='upper center',frameon=False,fontsize=10)
    save(fig,'figure4')
    return agg

def diagram():
    fig,ax=plt.subplots(figsize=(8.4,4.3));ax.set_xlim(0,10);ax.set_ylim(0,5);ax.axis('off')
    for x,title,subtitle in [(0.1,'Both groups represented','Query compatible with A and B'),(5.3,'Source group A omitted','Only the match to B remains')]:
        ax.text(x,4.7,title,fontweight='bold',fontsize=12)
        ax.text(x,4.28,subtitle,fontsize=10)
        ax.add_patch(FancyBboxPatch((x+.15,2.4),1.4,.65,boxstyle='round,pad=.12',facecolor='#E8EDF2',edgecolor='#556674'))
        ax.text(x+.85,2.72,'Query',ha='center',va='center',fontsize=11)
        for y,group in [(3.2,'Group A'),(1.6,'Group B')]:
            absent=x>5 and group=='Group A'
            ax.add_patch(FancyBboxPatch((x+2.7,y),1.4,.65,boxstyle='round,pad=.12',facecolor='white' if absent else '#D2E9E2',edgecolor='#8B9095',linestyle='--' if absent else '-'))
            ax.text(x+3.4,y+.32,'Omitted' if absent else group,ha='center',va='center',fontsize=10,color='#73787D' if absent else '#111111')
            if not absent:ax.add_patch(FancyArrowPatch((x+1.7,2.72),(x+2.53,y+.32),arrowstyle='->',mutation_scale=13,color='#007A5E',lw=1.7))
        ax.text(x,.75,'Ambiguous: no assignment' if x<5 else 'Neighbour match: diagnostic only',fontsize=11,fontweight='bold')
    ax.text(.1,.1,'Conceptual diagram of the observed failure mode; groups A and B are illustrative.',fontsize=9)
    save(fig,'figureS14')

def main():
    OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(exist_ok=True)
    parent=__import__('pathlib').Path(json.loads((RUN/'inputs/parent-run.json').read_text())['path'])
    for old,new in [('figure3','figureS11'),('figure4','figureS12')]:
        for ext in ['png','pdf','svg']:shutil.copyfile(parent/'outputs/figures'/f'{old}.{ext}',OUT/f'{new}.{ext}')
    result={'kraken':kraken(),'pterobacter':pterobacter(),'schematic':'S14 conceptual, not an additional experiment or an observed group assignment'}
    diagram()
    (RUN/'checks/v6-figure-data.json').write_text(json.dumps(result,indent=2)+'\n')
    for label,stems in [('manifest.json',[f'figure{i}' for i in range(1,5)]),('supplement-manifest.json',[f'figureS{i}' for i in range(1,15)])]:
        records=[]
        for stem in stems:
            for ext in ['png','pdf','svg']:
                p=OUT/f'{stem}.{ext}';records.append({'path':str(p.relative_to(RUN)),'sha256':sha(p),'bytes':p.stat().st_size})
        for p in sorted(DATA.glob('*')):records.append({'path':str(p.relative_to(RUN)),'sha256':sha(p),'bytes':p.stat().st_size})
        (OUT/label).write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
