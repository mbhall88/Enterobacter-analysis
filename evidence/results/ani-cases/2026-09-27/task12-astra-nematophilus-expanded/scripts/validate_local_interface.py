"""Pin assembly caller and exercise 6.25% divergence, deletion and reverse orientation."""
from pathlib import Path
import os,sys,random,json,subprocess,hashlib,shutil
r=Path(__file__).resolve().parents[1];p=r/'pinned';p.mkdir(exist_ok=True)
src=Path('${USER_HOME}/sw/snippy-ng/src/snippy_ng')
if not (p/'snippy_ng').exists():shutil.copytree(src,p/'snippy_ng',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
(p/'source-hashes.json').write_text(json.dumps({str(x.relative_to(p)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (p/'snippy_ng').rglob('*.py')},indent=2)+'\n')
d=r/'interface-test-v2';d.mkdir(exist_ok=True);rng=random.Random(20260927)
ref=''.join(rng.choices('ACGT',k=100000));expected=list(ref)
for i in sorted(rng.sample(range(len(ref)),6250)):expected[i]=rng.choice([b for b in 'ACGT' if b!=ref[i]])
query=''.join(expected[:40000]+expected[40100:]);query=query.translate(str.maketrans('ACGT','TGCA'))[::-1]
(d/'ref.fa').write_text('>ref\n'+ref+'\n');(d/'query.fa').write_text('>query\n'+query+'\n')
env=os.environ.copy();env['PATH']='${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin:'+env['PATH'];env['PYTHONPATH']=str(p)
for k,v in {'MUMMER_MIN_DELTA_IDENTITY':'0','MUMMER_MAX_DELTA_VARIANTS_PER_KB':'0','ASM_MAX_LOCAL_SNPS':'0'}.items():env['SNIPPY_NG_'+k]=v
cmd=['${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/snippy-ng','asm','--reference',str(d/'ref.fa'),'--assembly',str(d/'query.fa'),'--outdir',str(d/'out'),'--tmpdir',str(d/'tmp'),'--cpus','1','--ram','2','--aligner','nucmer','--no-cleanup']
(d/'command.json').write_text(json.dumps({'argv':cmd,'overrides':{k:v for k,v in env.items() if k.startswith('SNIPPY_NG_')},'source':'pinned/snippy_ng'},indent=2)+'\n')
with (d/'stdout').open('w') as o,(d/'stderr').open('w') as e:subprocess.run(cmd,env=env,stdout=o,stderr=e,check=True)
from Bio import SeqIO
seq=str(next(SeqIO.parse(d/'out/snippy.pseudo.fna','fasta')).seq).upper();assert len(seq)==len(ref)
callable=correct=0
for i,b in enumerate(seq):
 if 39980<=i<=40120:continue
 if b in 'ACGT':callable+=1;correct+=b==expected[i]
assert callable>98000 and correct/callable>.999, (callable,correct)
assert set(seq[40005:40095])<={'N','-'},set(seq[40005:40095])
(d/'validation.json').write_text(json.dumps({'status':'PASS','length':len(seq),'callable_excluding_indel_flanks':callable,'correct_callable':correct,'reverse_complement_and_deletion_checked':True,'designed_divergence':0.0625},indent=2)+'\n')
print('PASS',callable,correct)
