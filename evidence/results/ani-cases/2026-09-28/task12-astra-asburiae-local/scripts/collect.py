"""Collect completed evidence into tables and figures; author interpretation stays open."""
from pathlib import Path
import csv,json,hashlib,re,math
import numpy as np
from Bio import Phylo,SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from result_helpers import tips,key,splits,au_table
R=Path(__file__).resolve().parents[1];refs=['GCA_001521715.1','GCA_008364625.2']
def rows(p):return list(csv.DictReader((l for l in open(p) if not l.startswith('#')),delimiter='\t'))
def write(p,rs):
 with p.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
panel=rows(R/'inputs/panel.tsv');ids={x['accession'] for x in panel};groups={}
for x in panel:assert sha(x['source'])==x['sha256'];groups.setdefault(x['group'],set()).add(x['accession'])
for h in json.loads((R/'inputs/hypotheses.json').read_text()):groups[h['hypothesis']]=set(h['clade_members'])
groups['focal90']={x['accession'] for x in panel if x['focal']=='True'}
support=[];au=[];information=[];trees=[];edge_checks=0
for ref in refs:
 o=R/'local'/ref;t=o/'tests';status=json.loads((t/'validation.json').read_text());assert status['status']=='PASS'
 scfs=[]
 for seed in [20260928,20260929]:
  pref=t/('primary.scf-'+str(seed));stats={x['ID']:x for x in rows(Path(str(pref)+'.cf.stat'))};tr=Phylo.read(str(pref)+'.cf.branch','newick');d={}
  for k,n in splits(tr).items():
   branchid=n.name if n.name else str(int(n.confidence));assert branchid in stats;d[k]=stats[branchid]
  scfs.append(d)
 for label,path in [('primary',o/'tree/primary.treefile')]+([('RY',t/'RY.treefile')] if status['RY_run'] else []):
  tree=Phylo.read(path,'newick');assert tips(tree)==ids;sp=splits(tree);trees.append((ref,label,path))
  for name,ss in groups.items():
   k=key(ss,ids);n=sp.get(k);rec=dict(reference=ref,tree=label,group=name,tips=len(ss),monophyletic=n is not None,support_SHaLRT_UFBoot=(n.name or str(n.confidence)) if n is not None else '',branch_length=n.branch_length if n is not None else '',sCF_seed20260928='',sCF_seed20260929='')
   if label=='primary' and n is not None:
    for seed,d in zip([20260928,20260929],scfs):rec['sCF_seed'+str(seed)]=d[k]['sCF']
   support.append(rec)
 a=au_table((t/'topology/AU.iqtree').read_text());order=json.loads((t/'topology/candidate-order.json').read_text());assert len(a)==len({x['test_index'] for x in order})
 for x in order:au.append(dict(reference=ref,hypothesis=x['hypothesis'],test_index=x['test_index'],**a[x['test_index']]))
 for label,fasta in [('all93',o/'primary.core.fna'),('focal90',o/'focal90.core.fna')]:
  seqs={x.id:str(x.seq) for x in SeqIO.parse(fasta,'fasta')};arr=np.stack([np.frombuffer(s.encode(),dtype='S1') for s in seqs.values()]);assert set(seqs)==(ids if label=='all93' else groups['focal90']);counts=np.array([(arr==b.encode()).sum(axis=0) for b in 'ACGT']);states=(counts>0).sum(axis=0)
  information.append(dict(reference=ref,panel=label,taxa=len(seqs),columns=arr.shape[1],variable_sites=int((states>=2).sum()),parsimony_informative_sites=int(((counts>=2).sum(axis=0)>=2).sum()),invariant_sites=int((states==1).sum()),distinct_sequences=len(set(seqs.values()))))
write(R/'outputs/group-support.tsv',support);write(R/'outputs/topology-tests.tsv',au);write(R/'outputs/local-information.tsv',information)
# Figures show original distances; chinensis is the explicit display outgroup only.
COL={'asburiae_type_core':'#0072b2','dykesii_type_core':'#d55e00','shared_middle':'#009e73','peripheral_three':'#cc79a7','peripheral_six':'#946300','outside_boundary':'#555555','outside_boundary_partner':'#888888','chinensis_context':'#222222'}
LABEL={'asburiae_type_core':'A9','dykesii_type_core':'D3','shared_middle':'M69','peripheral_three':'P3','peripheral_six':'P6','outside_boundary':'boundary','outside_boundary_partner':'boundary partner','chinensis_context':'chinensis type'}
g={x['accession']:x['group'] for x in panel};(R/'figures').mkdir(exist_ok=False)
for kind in ['primary','RY']:
 selected=[x for x in trees if x[1]==kind]
 if not selected:continue
 fig,axes=plt.subplots(1,len(selected),figsize=(12*len(selected),22),squeeze=False)
 for ax,(ref,label,path) in zip(axes[0],selected):
  tr=Phylo.read(path,'newick');root='GCA_034008295.1';before={q:tr.distance(root,q) for q in ids};tr.root_with_outgroup(root);tr.ladderize()
  for q in ids:assert math.isclose(before[q],tr.distance(root,q),abs_tol=1e-10);edge_checks+=1
  leaves=tr.get_terminals();ys={n:float(i) for i,n in enumerate(leaves)};depth=tr.depths()
  for n in tr.get_nonterminals(order='postorder'):ys[n]=sum(ys[c] for c in n.clades)/len(n.clades)
  maxx=max(depth[n] for n in leaves);pad=maxx*.018
  for n in tr.find_clades():
   if n.clades:
    ax.plot([depth[n],depth[n]],[min(ys[c] for c in n.clades),max(ys[c] for c in n.clades)],color='#555555',lw=.6)
    for c in n.clades:ax.plot([depth[n],depth[c]],[ys[c],ys[c]],color=COL[g[c.name]] if c.is_terminal() else '#555555',lw=.7)
   else:
    text=n.name+'  '+LABEL[g[n.name]]+(' TYPE' if n.name in refs else '')
    ax.text(depth[n]+pad,ys[n],text,fontsize=7,color=COL[g[n.name]],va='center')
  model='nucleotide ML' if kind=='primary' else 'R/Y composition sensitivity'
  ax.set_title(('Asburiae' if ref==refs[0] else 'Dykesii')+' type reference\n'+model,fontsize=13);ax.set_ylim(-2,94);ax.invert_yaxis();ax.set_yticks([]);ax.set_xlim(0,maxx*1.9);ax.set_xlabel('Substitutions per site' if kind=='primary' else 'Substitutions per binary site');ax.spines[['top','right','left']].set_visible(False)
 fig.suptitle('Local 93-genome analysis: display groups are not assigned species',fontsize=16);fig.text(.5,.018,'Rooted on chinensis for display; tests use unrooted splits. P3/P6 mark prior placement disagreements.\nBranch support and two-seed site concordance: outputs/group-support.tsv. No recombination masking.',ha='center',fontsize=10);fig.tight_layout(rect=[0,.04,1,.97])
 for ext in ['png','pdf']:fig.savefig(R/'figures'/('local-'+kind+'.'+ext),dpi=180)
 plt.close(fig)
(R/'outputs/computational-validation.json').write_text(json.dumps(dict(status='PASS',references=2,input_hashes=93,support_rows=len(support),AU_rows=len(au),reroot_distance_checks=edge_checks,scientific_review='pending_agent_and_Michael',naming='unresolved',gubbins=False,dddh=False),indent=2)+'\n')
print('Computational collection complete; scientific interpretation and figure inspection pending.')
