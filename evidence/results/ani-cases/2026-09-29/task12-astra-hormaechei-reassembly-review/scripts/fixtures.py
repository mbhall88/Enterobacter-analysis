from pathlib import Path
import json,random,pysam
R=Path(__file__).resolve().parents[1];rng=random.Random(41962)
def read(p):
 with pysam.FastxFile(str(p)) as f:return {r.name:r.sequence for r in f}
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def err(s,rate):return ''.join(rng.choice([x for x in 'ACGT' if x!=c]) if rng.random()<rate else c for c in s)
models=read(R/'work/junction-models.fa');truth=[]
with (R/'work/span-controls.fa').open('w') as f:
 for model,s in models.items():
  for rate in [0,.001,.10]:
   for reverse in [False,True]:
    name=f'{model}_error{rate}_reverse{reverse}';seq=err(s[18500:53896],rate);seq=rc(seq) if reverse else seq
    f.write(f'>{name}\n{seq}\n');truth.append(dict(name=name,model=model,positive=True,substitution_rate=rate,reverse=reverse))
  for side,seq in [('left',s[:47000]),('right',s[25000:])]:
   name=model+'_'+side;f.write(f'>{name}\n{seq}\n');truth.append(dict(name=name,model=model,positive=False))
(R/'outputs/span-control-truth.json').write_text(json.dumps(truth,indent=2)+'\n')
# 50:50 repeat copies, ~40x per copy, 12kb reads. Mutated copy has five SNPs + 2bp deletion.
positions=[1000,8000,16320,25000,31000];variants=[]
with (R/'work/repeat-power.fq').open('w') as f:
 for copy,model in enumerate(['U1L-R-U2R','U2L-R-U1R']):
  s=models[model]
  if copy==1:
   t=list(s)
   for p in positions:
    old=t[20000+p];new=next(x for x in 'ACGT' if x!=old);t[20000+p]=new;variants.append(dict(repeat_pos1=p+1,reference=old,alternate=new))
   s=''.join(t);p=20000+20000;variants.append(dict(repeat_pos1=20000,reference=s[p-1],alternate='-'+s[p:p+2]));s=s[:p]+s[p+2:]
  for j,start in enumerate(range(0,len(s)-12000+1,300)):
   seq=err(s[start:start+12000],.001)
   if j%2:seq=rc(seq)
   f.write(f'@copy{copy}_{start}\n{seq}\n+\n'+('I'*len(seq))+'\n')
(R/'outputs/power-truth.json').write_text(json.dumps(dict(seed=41962,read_length=12000,step=300,substitution_rate=.001,variants=variants),indent=2)+'\n')
