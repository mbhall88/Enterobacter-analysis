"""Allele-specific, strand-aware checks; bases and indels reported separately."""
from pathlib import Path
from collections import Counter
import pysam,json,csv,statistics,sys
R=Path(__file__).resolve().parents[1];O=R.parent/'task12-opus55-hormaechei-reassembly'
def columns(bam,ref,ctg,start,end,bq=20):
 cache={};reference=ref.fetch(ctg)
 for c in bam.pileup(ctg,start,end,truncate=True,stepper='all',min_base_quality=0,min_mapping_quality=0,max_depth=100000,ignore_overlaps=False,ignore_orphans=False,flag_filter=0x704):
  base=Counter();indel=Counter();strand=Counter();mq0=0
  for p in c.pileups:
   a=p.alignment
   key=(a.query_name,a.flag,a.reference_start)
   if key not in cache:cache[key]=(a.query_sequence,a.query_qualities)
   seq,qual=cache[key]
   if p.is_refskip:continue
   if p.is_del:al='DEL'
   else:
    q=p.query_position
    if qual is not None and qual[q]<bq:continue
    al=seq[q].upper()
   base[al]+=1;strand[al+('(-)' if a.is_reverse else '(+)')]+=1;mq0+=(a.mapping_quality==0)
   if p.indel>0:indel['+'+seq[p.query_position+1:p.query_position+1+p.indel].upper()]+=1
   elif p.indel<0:indel['-'+reference[c.reference_pos+1:c.reference_pos+1-p.indel].upper()]+=1
  yield dict(contig=ctg,pos1=c.reference_pos+1,reference=reference[c.reference_pos],depth=sum(base.values()),mq0=mq0,bases=dict(base),indels=dict(indel),strands=dict(strand))
def scan(label,ctg,start,end,bq):
 bam=pysam.AlignmentFile(str(R/'work'/f'{label}.bam'));ref=pysam.FastaFile(str(R/'work'/f'{label.split("-sim")[0]}.fa'))
 rows=list(columns(bam,ref,ctg,start,end,bq));tested=[x for x in rows if x['depth']>=20];inter=[]
 for x in tested:
  alt={k:v for k,v in x['bases'].items() if k!=x['reference']};alt.update(x['indels'])
  for a,n in alt.items():
   if .2<=n/x['depth']<=.8:inter.append(dict(pos1=x['pos1'],allele=a,count=n,depth=x['depth'],frequency=n/x['depth']))
 out=dict(label=label,contig=ctg,start0=start,end0=end,base_quality_min=bq,positions_ge20=len(tested),positions_below20_or_uncovered=end-start-len(tested),depth_median=statistics.median([x['depth'] for x in rows]),mq0_fraction=sum(x['mq0'] for x in rows)/sum(x['depth'] for x in rows),intermediate_alleles=inter)
 (R/'outputs'/f'repeat-{label}-bq{bq}.json').write_text(json.dumps(out,indent=2)+'\n')
 with (R/'outputs'/f'repeat-{label}-bq{bq}.jsonl').open('w') as f:
  for x in rows:f.write(json.dumps(x)+'\n')
 bam.close();ref.close();return out
if __name__=='__main__':
 if '--sim' in sys.argv:
  for name,ctg,s,e in [('final-sim','chromosome_U1',4165525,4197921),('collapsed-sim','R_once',0,32396)]:print(json.dumps(scan(name,ctg,s,e,20)))
 else:
  for name,ctg,s,e in [('final','chromosome_U1',4165525,4197921),('collapsed','R_once',0,32396)]:
   for q in [0,20]:print(json.dumps(scan(name,ctg,s,e,q)))
  bam=pysam.AlignmentFile(str(R/'work/published.bam'));ref=pysam.FastaFile(str(R/'work/published.fa'));out=[]
  for line in (O/'outputs/published_hifi_disagree_sites.tsv').read_text().splitlines():
   c,p,*rest=line.split('\t')
   if not p.isdigit():continue
   out+=list(columns(bam,ref,c,int(p)-1,int(p),20))
  assert len(out)==47
  (R/'outputs/published-47-alleles.json').write_text(json.dumps(out,indent=2)+'\n')
  with (R/'outputs/published-47-alleles.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
