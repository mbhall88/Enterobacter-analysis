"""Check the original subreads BAM for long reads missing/changed in ENA FASTQ."""
from pathlib import Path
import pysam,json,hashlib,time
R=Path(__file__).resolve().parents[1]
D=Path('data/reads/hormaechei-type-SAMN16357575')
fq={}
with pysam.FastxFile(str(R/'work/subreads-ge20k.fq')) as f:
 for r in f:fq[(r.comment or r.name).split()[0]]=(len(r.sequence),hashlib.sha256(r.sequence.encode()).hexdigest())
n=long=0;missing=[];different=[];seen=set();start=time.time()
with pysam.AlignmentFile(str(D/'ccs/PFDA1_20200923_S64018_PL100159581_1_C01_bc1020_bc1020.bam'),'rb',check_sq=False) as b, (R/'work/additional-long-subreads.fq').open('w') as out:
 for r in b.fetch(until_eof=True):
  n+=1
  if r.query_length<20000:continue
  long+=1;name=r.query_name;seq=r.query_sequence;seen.add(name)
  if name not in fq:missing.append(name)
  elif fq[name]!=(len(seq),hashlib.sha256(seq.encode()).hexdigest()):different.append(name)
  else:continue
  qual=r.qual or ('I'*len(seq));out.write(f'@{name}\n{seq}\n+\n{qual}\n')
result=dict(bam_records=n,bam_ge20k=long,fastq_ge20k=len(fq),missing_in_fastq=missing,changed_sequence=different,fastq_only=sorted(set(fq)-seen),elapsed_seconds=time.time()-start)
(R/'outputs/bam-inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))
