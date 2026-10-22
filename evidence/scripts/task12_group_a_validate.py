#!/usr/bin/env python3
"""Independent checks of retained inputs, partitions and distance calculations."""
import argparse,csv,hashlib,json,itertools
from pathlib import Path
import numpy as np
from Bio import SeqIO
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score,silhouette_score
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from task12_group_a_populations import rows,write,sha,BASE,TYPES
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);out=ap.parse_args().output
checks=[]
def check(name,ok):
 checks.append(dict(check=name,passed=bool(ok)));assert ok,name
x=np.load(out/'inputs/matrices.npz');ids=list(x['ids']);n=len(ids);idx={q:i for i,q in enumerate(ids)};assign={r['accession']:int(r['operational_cluster']) for r in rows(out/'outputs/assignments.tsv')};labels=np.array([assign[q] for q in ids]);ut=np.triu_indices(n,1)
check('173 unique identifiers',len(set(ids))==n==173)
for line in (out/'inputs.sha256').read_text().splitlines():
 h,p=line.split('  ',1)
 check('unchanged input '+p,sha(p)==h)
check('complete pair ledger',len(rows(out/'inputs/pairs.tsv'))==n*(n-1)//2)
check('ANI symmetric, self100 and all AF80',np.array_equal(x['ani'],x['ani'].T) and np.all(np.diag(x['ani'])==100) and np.all(x['af']>=80))
for p in (out/'logs').glob('*.command.json'):
 r=json.loads(p.read_text());check('successful batch '+p.stem,r['exit_status']==0)
check('15 batch command records',len(list((out/'logs').glob('*.command.json')))==15)
for p in (out/'fastbaps').glob('*/clusters.tsv'):
 rr=rows(p);q=[r['accession'] for r in rr];expected=list(x['thin_ids']) if p.parent.name.startswith('thinned') else ids
 check('membership '+p.parent.name,len(q)==len(set(q)) and set(q)==set(expected))
# Raw pairwise-deletion counts against a direct byte comparison at selected pairs, independent of matrix multiplication.
sens=np.load(out/'coverage-sensitivity/pairwise-deletion.npz');check('95percent matrix identifiers',list(sens['ids'])==ids)
raw={r.id:np.frombuffer(str(r.seq).upper().encode(),dtype='S1') for r in SeqIO.parse(BASE/'core-genome/ledger-subset/alignment.fna','fasta') if r.id in set(ids)}
a=np.stack([raw[q] for q in ids]);valid=np.isin(a,np.frombuffer(b'ACGT',dtype='S1'));keep=valid.mean(0)>=.95
rng=np.random.default_rng(42);pairs=set(tuple(sorted(rng.choice(n,2,replace=False))) for _ in range(20));pairs.update(tuple(sorted((idx[q],idx[r]))) for q,r in itertools.combinations(TYPES,2))
for i,j in pairs:
 mask=keep&valid[i]&valid[j];nm=int(np.sum(a[i,mask]!=a[j,mask]));check(f'pairwise deletion counts {i},{j}',nm==sens['mismatches'][i,j] and int(mask.sum())==sens['callable'][i,j])
strict=np.load(out/'outputs/strict-core-distances.npy');rho=spearmanr(strict[ut],sens['distance'][ut]).statistic
extra=[]
for name in ['without_high_missing','intersection']:
 r=rows(out/f'coverage-sensitivity/{name}/clusters.tsv');pp={x['accession']:int(x['cluster']) for x in r};ii=sorted(pp);ari=adjusted_rand_score([assign[q] for q in ii],[pp[q] for q in ii]);extra.append(dict(dataset=name,n=len(ii),clusters=len(set(pp.values())),ARI_vs_primary=ari,types_separate=len({pp[q] for q in TYPES})==3))
write(out/'outputs/coverage-sensitivity.tsv',extra)
# Compare two distance-based hierarchies without selecting a biological rank or assuming k=3.
rr=[]
for k in [2,3,4,5,6,7,8,10,15,20]:
 aa=fcluster(linkage(squareform(strict),method='average'),k,criterion='maxclust');bb=fcluster(linkage(squareform(sens['distance']),method='average'),k,criterion='maxclust')
 rr.append(dict(max_clusters=k,strict_vs_pairwise_deletion_ARI=adjusted_rand_score(aa,bb),strict_clusters=len(set(aa)),pairwise_deletion_clusters=len(set(bb))))
write(out/'outputs/core-coverage-hierarchy-sensitivity.tsv',rr)
jac=np.load(out/'inputs/accessory-jaccard.npy');stats={'strict_vs_pairwise_deletion_spearman':rho,'pairwise_deletion_silhouette':silhouette_score(sens['distance'],labels,metric='precomputed'),'full_pangenome_Jaccard_silhouette':silhouette_score(jac,labels,metric='precomputed'),'additional_partitions':extra,'conditional_reference_alignment':'Not launched: population membership is assessed by existing ortholog alignment, separate global trees and near-full-core pairwise deletion. Remaining rank/hoffmannii residual issue is not a demonstrated reference-alignment defect.'}
(out/'outputs/coverage-summary.json').write_text(json.dumps(stats,indent=2)+'\n')
# Public cached metadata only; paired accessions taken explicitly from NCBI, not inferred by accession string replacement.
meta=Path('results/ani-cases/2026-10-02/task12-astra-sampling-audit/inputs/ncbi-enterobacter-current.jsonl');found={}
with meta.open() as f:
 for line in f:
  r=json.loads(line);matches=set([r['accession'],r.get('paired_accession','')])&set(ids)
  for q in matches:
   if q in found and r['accession']!=q:continue
   info=r.get('assembly_info',{});bs=info.get('biosample',{});attr={v['name']:v.get('value','') for v in bs.get('attributes',[])}
   found[q]=dict(accession=q,cluster=assign[q],metadata_accession=r['accession'],biosample=bs.get('accession',''),bioproject=info.get('bioproject_accession',''),host=attr.get('host',''),isolation_source=attr.get('isolation_source',''),geo_loc_name=attr.get('geo_loc_name',bs.get('geo_loc_name','')),strain=attr.get('strain',bs.get('strain','')))
write(out/'outputs/public-metadata.tsv',[found.get(q,dict(accession=q,cluster=assign[q],metadata_accession='',biosample='',bioproject='',host='',isolation_source='',geo_loc_name='',strain='')) for q in ids])
(out/'metadata-input.sha256').write_text(sha(meta)+'  '+str(meta)+'\n')
write(out/'outputs/metadata-coverage.tsv',[dict(cluster=c,n=sum(labels==c),matched=sum(q in found for q in ids if assign[q]==c),unique_biosamples=len({found[q]['biosample'] for q in found if assign[q]==c and found[q]['biosample']}),unique_bioprojects=len({found[q]['bioproject'] for q in found if assign[q]==c and found[q]['bioproject']})) for c in sorted(set(map(int,labels)))])
(out/'validation.json').write_text(json.dumps({'status':'PASS','checks':checks,'check_count':len(checks),'manual':'Figures reviewed; cluster rank not inferred; no new labels/assemblies; no formal phenotype association claimed.'},indent=2)+'\n');print(json.dumps(stats,indent=2));print('PASS',len(checks),'checks')
