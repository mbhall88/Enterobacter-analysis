from pathlib import Path
from collections import defaultdict,Counter
import pysam,json,re,hashlib
R=Path(__file__).resolve().parents[1];D=Path('data/reads/hormaechei-type-SAMN16357575')
def records(p):
 with pysam.FastxFile(str(p)) as f:return {x.name:(x.comment,x.sequence) for x in f}
# Tags: identical sequence inputs, only two header additions.
changes=[];inputs=[]
for p in sorted((D/'autocycler-hifi/assemblies').glob('*.fasta')):
 q=D/'autocycler-hifi-tagged/assemblies'/p.name
 assert q.exists();a,b=records(p),records(q);assert a.keys()==b.keys()
 for n in a:
  assert a[n][1]==b[n][1]
  if a[n][0]!=b[n][0]:
   assert b[n][0].replace('Autocycler_ignore','').strip()==(a[n][0] or '').strip()
   changes.append(dict(file=p.name,contig=n,before=a[n][0],after=b[n][0]))
 inputs.append(dict(file=p.name,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),tagged_sha256=hashlib.sha256(q.read_bytes()).hexdigest(),contigs=len(a)))
assert len(changes)==2
(R/'outputs/tag-check.json').write_text(json.dumps(dict(files=inputs,changes=changes),indent=2)+'\n')
# Real repeat-end allele counts, independently linked to >=1kb unique right flank.
assign=defaultdict(set)
for line in (R/'work/repeat-end-models.paf').open():
 f=line.split('\t')
 if int(f[7])<=52145 and int(f[8])>=53396:assign[f[0]].add(f[5].split('-')[-1])
b=pysam.AlignmentFile(str(R/'work/collapsed.bam'));groups=defaultdict(Counter);strand=defaultdict(Counter);perread=[]
for c in b.pileup('R_once',32145,32146,truncate=True,stepper='all',min_base_quality=20,max_depth=100000,ignore_overlaps=False):
 for p in c.pileups:
  if p.is_del or p.is_refskip:continue
  a=p.alignment;v='T8' if p.indel==0 else f'T{8+p.indel}'
  vals=assign.get(a.query_name,set());group=next(iter(vals)) if len(vals)==1 else ('ambiguous' if vals else 'unanchored')
  groups[group][v]+=1;strand[group][v+('-' if a.is_reverse else '+')]+=1
  perread.append(dict(read=a.query_name,right_flank=group,allele=v,mapq=a.mapping_quality))
(R/'outputs/homopolymer-flank-check.json').write_text(json.dumps(dict(groups={k:dict(v) for k,v in groups.items()},strands={k:dict(v) for k,v in strand.items()},reads=perread),indent=2)+'\n')
# Is the independent2016 scaffold bridge gap-free? Check sequence within actual alignment.
draft=records(R/'type-census/raw/GCA_001875655.1_genomic.fna.gz');out=[]
for line in (R/'work/GCA_001875655.1-models.paf').open():
 f=line.split('\t')
 if int(f[7])<=19000 and int(f[8])>=53396:
  s=draft[f[0]][1][int(f[2]):int(f[3])];out.append(dict(scaffold=f[0],query_start0=int(f[2]),query_end0=int(f[3]),model=f[5],matches=int(f[9]),alignment_columns=int(f[10]),N_bases=s.upper().count('N'),ambiguous_bases=sum(c.upper() not in 'ACGT' for c in s)))
(R/'outputs/independent-draft-bridge.json').write_text(json.dumps(out,indent=2)+'\n')
print('tag and independent draft checks pass');print(json.dumps({k:dict(v) for k,v in groups.items()}));print(json.dumps(out))
