"""Build validated coordinate-preserving local alignments; keep all invariant core sites."""
import argparse,csv,json,os,subprocess,time,hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from Bio import SeqIO
from mask_pseudo import masked_sequences
r=Path(__file__).resolve().parents[1];root=Path('${PROJECT_ROOT}')
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--threads',type=int,default=4);a=p.parse_args()
rows=list(csv.DictReader((r/'inputs/panel.tsv').open(),delimiter='\t'));byid={x['accession']:x for x in rows};assert a.reference in byid
out=r/'local'/a.reference;out.mkdir(parents=True,exist_ok=False)
assert json.loads((r/'outputs/preflight.json').read_text())['status']=='PASS'
records=list(SeqIO.parse(byid[a.reference]['source'],'fasta'));ref={x.id:str(x.seq).upper() for x in records}
# Keep contig boundaries explicit with uncallable padding and a coordinate ledger.
padding=10000;offset=0;coords=[]
for i,(name,seq) in enumerate(ref.items()):coords.append(dict(contig=name,start0=offset,end0=offset+len(seq),padding_after=padding if i<len(ref)-1 else 0));offset+=len(seq)+(padding if i<len(ref)-1 else 0)
with (out/'coordinates.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(coords[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(coords)
env=os.environ.copy();env['PATH']='${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin:'+env['PATH'];env['PYTHONPATH']=str(r/'pinned');env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
for k,v in {'MUMMER_MIN_DELTA_IDENTITY':'0','MUMMER_MAX_DELTA_VARIANTS_PER_KB':'0','ASM_MAX_LOCAL_SNPS':'0'}.items():env['SNIPPY_NG_'+k]=v
selected=rows;assert len(selected)==93
for x in selected:assert hashlib.sha256(Path(x['source']).read_bytes()).hexdigest()==x['sha256']
(out/'overrides.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith('SNIPPY_NG_')},indent=2)+'\n')
def work(x):
 q=x['accession']
 if q==a.reference:return q,ref,None
 d=out/'samples'/q;d.parent.mkdir(exist_ok=True)
 cmd=['${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/snippy-ng','asm','--reference',byid[a.reference]['source'],'--assembly',x['source'],'--outdir',str(d),'--tmpdir',str(out/'tmp'/q),'--cpus','1','--ram','2','--aligner','nucmer','--no-cleanup']
 (out/'logs').mkdir(exist_ok=True);before=time.monotonic()
 with (out/'logs'/(q+'.stdout')).open('w') as o,(out/'logs'/(q+'.stderr')).open('w') as e:
  status=subprocess.run(['/usr/bin/time','-v','-o',str(out/'logs'/(q+'.time'))]+cmd,stdout=o,stderr=e,env=env).returncode
 rec=dict(accession=q,argv=cmd,seconds=time.monotonic()-before,exit_code=status,slurm_job=os.getenv('SLURM_JOB_ID'));(out/'logs'/(q+'.command.json')).write_text(json.dumps(rec,indent=2)+'\n');assert status==0,q
 seqs=masked_sequences(d);assert list(seqs)==list(ref) and all(len(seqs[k])==len(ref[k]) for k in ref)
 print(q,round(rec['seconds'],1),flush=True);return q,seqs,rec
aligned={};stats=[]
with ThreadPoolExecutor(max_workers=a.threads) as pool:
 for q,seqs,rec in pool.map(work,selected):
  s=('N'*padding).join(seqs[k] for k in ref);assert len(s)==offset
  aligned[q]=np.frombuffer(s.encode(),dtype='S1').copy();stats.append(dict(accession=q,callable_bases=sum(s.count(b) for b in 'ACGT'),reference_genomic_bases=sum(map(len,ref.values())),coordinate_length=offset,focal=byid[q]['focal']))
ids=sorted(aligned);matrix=np.stack([aligned[q] for q in ids]);isbase=np.isin(matrix,np.array(list('ACGT'),dtype='S1'))
primary=ids;assert len(primary)==93
ix=[ids.index(q) for q in primary];core=isbase[ix].all(axis=0)
assert core.sum()>1000000, int(core.sum())
np.save(out/'strict-core-mask.npy',core)
def fasta(path,names,positions=None,mask=False):
 with path.open('w') as f:
  for q in names:
   seq=aligned[q].copy()
   if mask:seq[~core]=b'N'
   if positions is not None:seq=seq[positions]
   f.write('>'+q+'\n'+seq.tobytes().decode()+'\n')
fasta(out/'all-projected.fna',ids)
fasta(out/'primary.full-coordinate.fna',primary,mask=True)
fasta(out/'primary.core.fna',primary,positions=core)
# Quantify extra context's impact on retained sites without silently changing the panel.
focal=[q for q in ids if byid[q]['focal']=='True'];assert len(focal)==90
focal_core=isbase[[ids.index(q) for q in focal]].all(axis=0)
np.save(out/'focal-only-core-mask.npy',focal_core)
fasta(out/'focal90.core.fna',focal,positions=focal_core)
with (out/'callable.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(stats[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(stats)
sub=matrix[ix][:,core];variable=(sub!=sub[0]).any(axis=0)
(out/'validation.json').write_text(json.dumps({'status':'PASS','reference':a.reference,'primary_genomes':93,'projected_genomes':93,'focal90_core_bases':int(focal_core.sum()),'coordinate_length':offset,'reference_bp':sum(map(len,ref.values())),'strict_core_bases':int(core.sum()),'variable_core_sites':int(variable.sum()),'invariant_core_sites':int((~variable).sum()),'padding':padding,'gubbins':'not_run','indel_mask_padding':10,'identity_or_snp_density_filter':'disabled; divergence must not preselect taxonomy','rows':ids,'commands':'logs/*.command.json'},indent=2)+'\n')
print(json.loads((out/'validation.json').read_text()),flush=True)
