from pathlib import Path
from collections import defaultdict
import re,csv,json,pysam
R=Path(__file__).resolve().parents[1]
def analyse(path):
 groups=defaultdict(list);hits=[]
 for line in path.open():
  f=line.rstrip().split('\t');q,model=f[0],f[5];ts,te=int(f[7]),int(f[8]);ql=int(f[1]);cg=next(x[5:] for x in f[12:] if x.startswith('cg:Z:'));t=ts;left=right=0;gap=0
  for nn,op in re.findall(r'(\d+)([MIDNSHP=X])',cg):
   n=int(nn)
   if op in 'M=X':left+=max(0,min(t+n,20000)-t);right+=max(0,t+n-max(t,52396))
   if op in 'DN' and t<52396 and t+n>20000:gap=max(gap,n)
   if op in 'MDN=X':t+=n
  assert t==te
  span=left>=1000 and right>=1000 and gap<=1000
  extent=ts<=19000 and te>=53396
  reach=max(min(te,52396)-20000 if ts<=19000 else 0,52396-max(ts,20000) if te>=53396 else 0,0)
  row=dict(read=q,query_length=ql,model=model,strand=f[4],query_start=int(f[2]),query_end=int(f[3]),target_start=ts,target_end=te,left_aligned=left,right_aligned=right,max_repeat_deletion=gap,spans=span,extent_spans=extent,anchored_repeat_reach=reach)
  hits.append(row);groups[q,model,f[4]].append(row)
 envelopes=[dict(read=k[0],model=k[1],strand=k[2],hits=len(v)) for k,v in groups.items() if min(x['target_start'] for x in v)<=19000 and max(x['target_end'] for x in v)>=53396]
 return dict(paf=str(path),alignment_count=len(hits),reads_with_alignments=len({x['read'] for x in hits}),spanning_hits=[x for x in hits if x['spans']],extent_spanning_hits=[x for x in hits if x['extent_spans']],envelope_candidates=envelopes,best_reach=sorted(hits,key=lambda x:x['anchored_repeat_reach'],reverse=True)[:12])
for label in ['hifi-allchains','subreads-allchains','controls']:
 out=analyse(R/'work'/f'span-{label}.paf')
 if label!='controls':
  fp='hifi-ge20k.fq' if label.startswith('hifi') else 'subreads-ge20k.fq'
  with pysam.FastxFile(str(R/'work'/fp)) as f:records=[(r.name,r.comment,len(r.sequence)) for r in f]
  out['input_reads_ge20k']=len(records);out['input_reads_ge30k']=sum(x[2]>=30000 for x in records)
  out['input_zmws_ge20k']=len({(x[1] or x[0]).rsplit('/',1)[0] for x in records})
 else:
  truth=json.loads((R/'outputs/span-control-truth.json').read_text());spans={x['read'] for x in out['spanning_hits']}
  out['positive_controls']=sum(x['positive'] for x in truth);out['negative_controls']=sum(not x['positive'] for x in truth)
  out['missed_positives']=[x for x in truth if x['positive'] and x['name'] not in spans];out['false_positives']=[x for x in truth if not x['positive'] and x['name'] in spans]
 (R/'outputs'/f'span-{label}.json').write_text(json.dumps(out,indent=2)+'\n')
 print(label,json.dumps({k:v for k,v in out.items() if k!='best_reach'}))
for a in ['GCA_001875655.1','GCA_000213995.1']:
 p=R/'work'/f'{a}-models.paf'
 if p.exists():
  out=analyse(p);(R/'outputs'/f'{a}-models.json').write_text(json.dumps(out,indent=2)+'\n')
