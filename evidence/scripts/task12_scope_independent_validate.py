#!/usr/bin/env python3
"""Check interval unions with an independent endpoint sweep and deletion bounds."""
import collections,json
import task12_scope_s1_s2 as s
def sweep(v):
 e=collections.Counter()
 for a,b in v:e[min(a,b)]+=1;e[max(a,b)+1]-=1
 depth=total=0;prev=None
 for pos,change in sorted(e.items()):
  if depth:total+=pos-prev
  depth+=change;prev=pos
 return total
n=0
for j in s.rows(s.RUN/'alignment-jobs.tsv'):
 target=s.fasta(j['target']);cats=collections.defaultdict(list);p=s.RUN/'alignments'/(j['kind']+'_'+j['accession'])/'neighbour.coords'
 with open(p) as f:
  for line in f:
   v=line.rstrip().split('\t')
   if float(v[6])>=95:cats[(v[10],s.label(target[v[9]][0]))].append((int(v[2]),int(v[3])))
 for v in cats.values():assert sweep(v)==s.bp(v);n+=1
bounds=s.rows(s.RUN/'s2-pair-selection-census.tsv');a={r['accession']:r for r in s.rows(s.RUN/'s2-anchors.tsv')}
for r in bounds:
 lo,hi=s.deletion_bounds(float(r['original_skani_type_af']),int(a[r['type_accession']]['draft_bp']),int(r['removed_bp']))
 assert abs(lo-float(r['fixed_alignment_type_af_lower_pct']))<1e-12 and abs(hi-float(r['fixed_alignment_type_af_upper_pct']))<1e-12
(s.RUN/'independent-validation.json').write_text(json.dumps(dict(event_sweep_union_groups=n,all_match=True,deletion_bounds_checked=len(bounds)),indent=2));print(n)
