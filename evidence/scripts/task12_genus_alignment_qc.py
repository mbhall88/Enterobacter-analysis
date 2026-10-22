#!/usr/bin/env python3
import argparse,csv,json
from collections import Counter
from pathlib import Path
from task12_prepare_genus_tree import fasta,writefa,writetsv,AA
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run
stats=[];dups=[]
for label in ['primary-full80','canonical-masked']:
 seqs=dict(fasta(r/'inputs'/f'{label}.faa'));var=pi=empty=0
 for col in zip(*seqs.values()):
  c=Counter(x for x in col if x in AA);var+=len(c)>1;pi+=sum(n>=2 for n in c.values())>=2;empty+=not c
 groups={}
 for n,s in seqs.items():groups.setdefault(s,[]).append(n)
 for members in groups.values():
  if len(members)>1:dups.append(dict(alignment=label,count=len(members),tips=';'.join(sorted(members))))
 stats.append(dict(alignment=label,tips=len(seqs),columns=len(next(iter(seqs.values()))),variable=var,parsimony_informative=pi,empty=empty,constant=len(next(iter(seqs.values())))-var,distinct_sequences=len(groups)))
assert stats[0]['parsimony_informative']==6644 and stats[1]['parsimony_informative']==206
writetsv(r/'tables/alignment-statistics.tsv',stats);writetsv(r/'tables/identical-alignments.tsv',dups)
full=dict(fasta(r/'inputs/primary-full80.faa'));assert 'GCA_023425455.1' in full;del full['GCA_023425455.1']
writefa(r/'inputs/primary-full80-no-composition-outlier.faa',full)
(r/'inputs/composition-sensitivity.json').write_text(json.dumps(dict(excluded_reference='GCA_023425455.1',reason='Only primary-alignment sequence failing IQ-TREE composition chi-square; p=0.0002;33.63%gap/ambiguity. Diagnostic sensitivity only, not a biological exclusion.',kept_columns='identical to primary-full80; no re-filtering',tips=len(full)),indent=2)+'\n')
print(stats)
