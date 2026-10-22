#!/usr/bin/env python3
"""Audit existing Group A data and summarize label-independent population partitions."""
import argparse, csv, hashlib, itertools, json, platform, sys
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform, pdist
from Bio import SeqIO, Phylo

TYPES={'GCA_001729745.1':'hoffmannii','GCA_040096145.1':'intestinihominis','GCA_052922205.1':'nototheniae'}
BASE=Path('results/conflict-resolution/2026-09-23/task17-opus55-resolve/C3-trio')
MEM=Path('results/ani-cases/2026-09-28/task12-astra-ab-review/outputs-v2/genomes.tsv')
PAIRS=Path('results/ani-cases/2026-09-23/task12-astra-evidence/skani/pairs-complete.tsv')
CFG=Path('results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/tree-comparison.json')
def rows(p):
 with open(p) as f: return list(csv.DictReader(f,delimiter='\t'))
def write(p,rr):
 rr=list(rr)
 if not rr: return
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def canonical(labels,ids):
 groups={c:sorted(q for q,l in zip(ids,labels) if l==c) for c in set(labels)}
 order={c:i+1 for i,c in enumerate(sorted(groups,key=lambda c:groups[c][0]))}
 return np.array([order[c] for c in labels])
def fasta(path,ids,m):
 with open(path,'w') as f:
  for q,x in zip(ids,m):f.write('>'+q+'\n'+x.tobytes().decode()+'\n')
def strict_core(a):
 valid=np.isin(a,np.frombuffer(b'ACGT',dtype='S1')); keep=np.all(valid,axis=0)
 core=a[:,keep];var=np.any(core!=core[0],axis=0)
 return valid,keep,var

def prepare(out):
 out.mkdir(exist_ok=True); inp=out/'inputs';inp.mkdir(exist_ok=False)
 members={r['accession']:r for r in rows(MEM) if r['union_group']=='A'};ids=sorted(members); n=len(ids)
 assert n==173 and set(TYPES)<=set(ids)
 subset={r['tip_id']:r for r in rows(BASE/'subset/subset.tsv')}; assert set(ids)<=set(subset) and len(subset)==178
 aln=BASE/'core-genome/ledger-subset/alignment.fna';seq={r.id:np.frombuffer(str(r.seq).upper().encode(),dtype='S1') for r in SeqIO.parse(aln,'fasta')}
 assert set(seq)==set(subset) and len({len(s) for s in seq.values()})==1
 a=np.stack([seq[q] for q in ids]);del seq
 valid,keep,var=strict_core(a);ix=np.where(keep)[0][var];snps=a[:,ix]
 assert snps.shape[1]>0
 ledger=rows(BASE/'core-genome/ledger-subset/ledger.tsv');ledger=[r for r in ledger if r['in_concatenation']=='True']
 gene=np.full(a.shape[1],-1,dtype=np.int32)
 for i,r in enumerate(ledger):
  start,end=int(r['start'])-1,int(r['end']);assert np.all(gene[start:end]==-1);gene[start:end]=i
 assert (gene>=0).all();sgene=gene[ix]
 ani=np.full((n,n),np.nan);np.fill_diagonal(ani,100);af=np.full((n,n),np.nan);np.fill_diagonal(af,100);idx={q:i for i,q in enumerate(ids)};pairrows=[]
 with PAIRS.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['a'] in idx and r['b'] in idx:
    i,j=idx[r['a']],idx[r['b']];assert np.isnan(ani[i,j]) and r['status']=='measured'
    ani[i,j]=ani[j,i]=float(r['ani']);af[i,j]=af[j,i]=min(float(r['af_a']),float(r['af_b']));pairrows.append(r)
 assert len(pairrows)==n*(n-1)//2 and np.isfinite(ani).all()
 # Complete linkage limits maximum within-group ANI divergence; representative thinning is not isolate deduplication.
 z=linkage(squareform(100-ani,checks=True),method='complete');thinlabels=fcluster(z,t=.05+1e-10,criterion='distance');selected=set(TYPES)
 for c in set(thinlabels):
  ss=[q for q,l in zip(ids,thinlabels) if l==c]
  if not selected.intersection(ss):selected.add(min(ss))
 thinidx=[i for i,q in enumerate(ids) if q in selected]
 fasta(inp/'primary.fasta',ids,snps);fasta(inp/'thinned.fasta',[ids[i] for i in thinidx],snps[thinidx])
 design=[{'dataset':'primary','prior':p,'alignment':str(inp/'primary.fasta')} for p in ['optimised.symmetric','symmetric','baps','optimised.baps']]
 design.append({'dataset':'thinned','prior':'optimised.symmetric','alignment':str(inp/'thinned.fasta')})
 rng=np.random.default_rng(20261002)
 for rep in range(5):
  selected_genes=rng.permutation(len(ledger))[:len(ledger)//2];chosen=np.isin(sgene,selected_genes)
  for side,mask in [('a',chosen),('b',~chosen)]:
   name=f'half{rep+1}{side}';fasta(inp/f'{name}.fasta',ids,snps[:,mask]);design.append({'dataset':name,'prior':'optimised.symmetric','alignment':str(inp/f'{name}.fasta')})
 np.savez_compressed(inp/'matrices.npz',ani=ani,af=af,snps=snps,site_gene=sgene,ids=np.array(ids),strict_core_columns=int(keep.sum()),thin_ids=np.array(sorted(selected)))
 write(inp/'datasets.tsv',design);write(inp/'pairs.tsv',pairrows)
 write(inp/'membership.tsv',[dict(accession=q,type_name=TYPES.get(q,''),archived_label=members[q]['archived_label'],gtdb_species=subset[q]['gtdb_species'],in_panel=subset[q]['in_panel'],in_all_three_minimal_type_clades=members[q]['in_all_three'],callable_bases=int(valid[i].sum()),alignment_bases=a.shape[1],missing_fraction=float(1-valid[i].mean()),in_thinned=q in selected) for i,q in enumerate(ids)])
 retained_per_gene=np.bincount(gene[keep],minlength=len(ledger)); variable_per_gene=np.bincount(sgene,minlength=len(ledger))
 write(inp/'locus-retention.tsv',[dict(locus=r['locus'],original_columns=int(r['end'])-int(r['start'])+1,strict_core_columns=int(retained_per_gene[i]),variable_columns=int(variable_per_gene[i])) for i,r in enumerate(ledger)])
 # Complementary gene-content distances: inherited whole C3 pangenome, restricted to A.
 gp=BASE/'core-genome/panaroo/gene_presence_absence.Rtab'
 with gp.open() as f:
  reader=csv.reader(f,delimiter='\t');head=next(reader);positions=[head.index(q) for q in ids];content=np.array([[int(r[k])>0 for k in positions] for r in reader],dtype=bool).T
 jac=squareform(pdist(content,metric='jaccard'));np.save(inp/'accessory-jaccard.npy',jac)
 summary={'n':n,'intersection_minimal_type_clades':sum(r['in_all_three']=='True' for r in members.values()),'non_A_excluded':sorted(set(subset)-set(ids)),'alignment_columns':a.shape[1],'strict_core_columns':int(keep.sum()),'strict_core_fraction':float(keep.mean()),'variable_columns':snps.shape[1],'loci':len(ledger),'loci_with_variable_columns':len(set(sgene)),'thinned_n':len(selected),'seed':20261002,'datasets':len(design),'all_pairs_af_ge80':bool((af>=80).all()),'original_pan_genes':content.shape[1],'variable_pan_genes':int(((content.sum(0)>0)&(content.sum(0)<n)).sum())}
 (out/'preparation.json').write_text(json.dumps(summary,indent=2)+'\n')
 sources=[MEM,PAIRS,CFG,aln,BASE/'subset/subset.tsv',BASE/'core-genome/ledger-subset/ledger.tsv',gp,Path(__file__)]
 (out/'inputs.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in sources))
 print(json.dumps(summary,indent=2))

def describe(x):
 return dict(n=len(x),min=float(np.min(x)) if len(x) else '',median=float(np.median(x)) if len(x) else '',max=float(np.max(x)) if len(x) else '')
def tree_metrics(path,groups,ids):
 t=Phylo.read(path,'newick');tips={x.name:x for x in t.get_terminals()};universe=set(tips)&set(ids);rr=[]
 for c in sorted(set(groups.values())):
  ss={q for q in universe if groups[q]==c}
  if not ss:continue
  node=t.common_ancestor([tips[q] for q in sorted(ss)]);desc={x.name for x in node.get_terminals()};other=(desc & universe)-ss
  rr.append(dict(cluster=c,represented=len(ss),full_cluster_size=sum(v==c for v in groups.values()),mrca_all_tips=len(desc),other_A_in_mrca=len(other),pure_in_A=not other,monophyletic_all_tips=desc==ss,support=node.confidence if node.confidence is not None else node.name or '',branch_length=node.branch_length or 0))
 return rr

def summarize(out):
 from sklearn.metrics import adjusted_rand_score,silhouette_score
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 d=np.load(out/'inputs/matrices.npz');ids=list(d['ids']);n=len(ids);idx={q:i for i,q in enumerate(ids)};ani=d['ani'];snps=d['snps'];dist=100-ani
 res=out/'outputs';res.mkdir(exist_ok=False);fig=out/'figures';fig.mkdir(exist_ok=False)
 parts={}
 for p in sorted((out/'fastbaps').glob('*/clusters.tsv')):
  rr=rows(p);labs={r['accession']:int(r['cluster']) for r in rr};assert len(labs)==len(rr) and set(labs)<=set(ids);parts[p.parent.name]=labs
 assert len(parts)==15,len(parts)
 primary=parts['primary--optimised.symmetric'];assert set(primary)==set(ids);labels=canonical([primary[q] for q in ids],ids);primary=dict(zip(ids,map(int,labels)))
 cfg=json.loads(CFG.read_text());treerows=[]
 for method,v in cfg['comparison_sets']['ingroup-1076']['rooting']['per_tree'].items():
  p=v['rootings']['both_set12_outgroups']['tree'];treerows.extend(dict(method=method,scope='full_global',**r) for r in tree_metrics(p,primary,ids))
 for method,p in [('core_all_positions',BASE/'core-genome/topology-panel-codon3-fixedparams/au.treefile'),('core_first_second',BASE/'core-genome/topology-panel-pos12-fixedparams/au.treefile')]:
  treerows.extend(dict(method=method,scope='53_genome_panel_restriction',**r) for r in tree_metrics(p,primary,ids))
 write(res/'tree-cluster-metrics.tsv',treerows)
 members={r['accession']:r for r in rows(out/'inputs/membership.tsv')};write(res/'assignments.tsv',[dict(**members[q],operational_cluster=primary[q]) for q in ids])
 sensitivity=[];co=[]
 for name,part in parts.items():
  ii=[q for q in ids if q in part];ll=[part[q] for q in ii];sizes=[ll.count(c) for c in set(ll)];ss={q:part[q] for q in TYPES};same=ss[list(TYPES)[0]]==ss[list(TYPES)[1]]==ss[list(TYPES)[2]]
  sensitivity.append(dict(analysis=name,n=len(ii),clusters=len(set(ll)),singletons=sum(x==1 for x in sizes),largest=max(sizes),adjusted_rand_vs_primary=adjusted_rand_score([primary[q] for q in ii],ll),types_same_cluster=same,**{TYPES[q]+'_cluster_size':ll.count(part[q]) for q in TYPES}))
  if name.startswith('half'):co.append(np.equal.outer([part[q] for q in ids],[part[q] for q in ids]))
 co=np.mean(co,axis=0);np.save(res/'gene-half-coassignment.npy',co);write(res/'sensitivity.tsv',sensitivity)
 type_pairs=[]
 for a,b in itertools.combinations(TYPES,2):
  type_pairs.append(dict(a=TYPES[a],b=TYPES[b],ani=ani[idx[a],idx[b]],same_primary=primary[a]==primary[b],half_gene_coassignment=co[idx[a],idx[b]],same_priors=';'.join(name for name,p in parts.items() if name.startswith('primary--') and p[a]==p[b])))
 write(res/'type-pairs.tsv',type_pairs)
 between=[]
 for ca in sorted(set(labels)):
  for cb in sorted(set(labels)):
   if ca>cb:continue
   ia=np.where(labels==ca)[0];ib=np.where(labels==cb)[0];v=ani[np.ix_(ia,ib)];x=v[np.triu_indices(len(ia),1)] if ca==cb else v.ravel()
   between.append(dict(cluster_a=int(ca),cluster_b=int(cb),kind='within' if ca==cb else 'between',**describe(x)))
 write(res/'cluster-ANI.tsv',between)
 nearest=[]
 for i,q in enumerate(ids):
  outsiders=[j for j in range(n) if labels[j]!=labels[i]]
  if not outsiders:continue
  best=max(ani[i,j] for j in outsiders);hits=[ids[j] for j in outsiders if ani[i,j]==best]
  within=[ani[i,j] for j in range(n) if j!=i and labels[j]==labels[i]]
  nearest.append(dict(accession=q,cluster=primary[q],nearest_outside_ani=best,tied_nearest_outside=';'.join(hits),within_min=min(within) if within else '',within_max=max(within) if within else '',outside_exceeds_farthest_within=bool(within and best>min(within))))
 write(res/'nearest-outside.tsv',nearest)
 thresholds=[]
 for method in ['average','complete','single']:
  z=linkage(squareform(dist),method=method)
  for cutoff in [98,98.5,99,99.2,99.5,99.8,99.9,99.95]:
   l=fcluster(z,100-cutoff+1e-10,criterion='distance');sizes=[int(sum(l==c)) for c in set(l)]
   thresholds.append(dict(method=method,cutoff=cutoff,clusters=len(sizes),largest=max(sizes),**{TYPES[q]+'_size':int(sum(l==l[idx[q]])) for q in TYPES},ARI_vs_fastbaps=adjusted_rand_score(labels,l)))
 write(res/'ANI-threshold-sensitivity.tsv',thresholds)
 core=squareform(pdist(snps.view('u1'),metric='hamming'))*snps.shape[1]/int(d['strict_core_columns']);np.save(res/'strict-core-distances.npy',core)
 write(res/'core-distances.tsv',[dict(a=ids[i],b=ids[j],core_mismatches=round(core[i,j]*int(d['strict_core_columns'])),callable_core_bases=int(d['strict_core_columns']),core_divergence=core[i,j],skani_ani=ani[i,j],same_cluster=labels[i]==labels[j]) for i,j in itertools.combinations(range(n),2)])
 overview=[]
 for c in sorted(set(labels)):
  ii=np.where(labels==c)[0];ss=[ids[i] for i in ii];off=co[np.ix_(ii,ii)][np.triu_indices(len(ii),1)]
  overview.append(dict(cluster=int(c),n=len(ii),types=';'.join(TYPES[q] for q in ss if q in TYPES),half_gene_within_mean=float(off.mean()) if len(off) else '',half_gene_within_min=float(off.min()) if len(off) else '',source_labels=';'.join(sorted(set(members[q]['archived_label'] for q in ss)))))
 write(res/'clusters.tsv',overview)
 order=leaves_list(linkage(squareform(dist),method='average'));fig1,ax=plt.subplots(1,2,figsize=(14,6))
 im=ax[0].imshow(ani[np.ix_(order,order)],vmin=98.5,vmax=100,cmap='viridis');fig1.colorbar(im,ax=ax[0],label='Chromosome skani ANI (%)')
 im=ax[1].imshow(co[np.ix_(order,order)],vmin=0,vmax=1,cmap='magma');fig1.colorbar(im,ax=ax[1],label='Same cluster in ten half-gene analyses')
 tp=[int(np.where(order==idx[q])[0][0]) for q in TYPES]
 for axis,title in zip(ax,['Existing ANI: all 173 Group A genomes','Population stability: random gene halves']):
  axis.set_xticks(tp,list(TYPES.values()),rotation=45,ha='right');axis.set_yticks(tp,list(TYPES.values()));axis.set_title(title)
 fig1.tight_layout();fig1.savefig(fig/'ani-and-stability.png',dpi=180);plt.close(fig1)
 fig2,ax=plt.subplots(1,2,figsize=(13,5));ut=np.triu_indices(n,1);same=np.equal.outer(labels,labels)[ut]
 for mask,label,color in [(same,'Within fastbaps populations','#0072B2'),(~same,'Between fastbaps populations','#D55E00')]:
  ax[0].hist(ani[ut][mask],bins=np.linspace(98.5,100,61),alpha=.55,label=label,color=color)
 ax[0].set(xlabel='Chromosome skani ANI (%)',ylabel='Pair count',title='Do inferred populations have an ANI gap?');ax[0].legend(fontsize=8)
 for k,(q,name) in enumerate(TYPES.items()):
  x=ani[idx[q],:];o=np.argsort(x);ax[1].plot(np.arange(n-1),np.sort(np.delete(x,idx[q])),label=name)
 ax[1].set(xlabel='Other genomes, ordered separately for each type',ylabel='Chromosome skani ANI (%)',title='Similarity to each Group A type');ax[1].legend()
 fig2.tight_layout();fig2.savefig(fig/'population-ani-distributions.png',dpi=180);plt.close(fig2)
 summary={'n':n,'primary_clusters':len(set(labels)),'primary_sizes':sorted([int(sum(labels==c)) for c in set(labels)],reverse=True),'silhouette_ANI':float(silhouette_score(dist,labels,metric='precomputed')) if 1<len(set(labels))<n else None,'silhouette_core':float(silhouette_score(core,labels,metric='precomputed')) if 1<len(set(labels))<n else None,'type_pairs':type_pairs,'sensitivity':sensitivity,'strict_core_bases':int(d['strict_core_columns']),'no_species_or_subspecies_assignments_changed':True}
 (res/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['prepare','summarize']);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();{'prepare':prepare,'summarize':summarize}[a.stage](a.output)
if __name__=='__main__':main()
