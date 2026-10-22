#!/usr/bin/env python3
"""Describe profiles by build-accession evidence and study type-anchor role."""
import csv,collections,statistics
from pathlib import Path
R=Path(__file__).resolve().parent

def read(path):
 with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))
def write(path,rows):
 with open(path,'w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
account={r['accession']:r for r in read(R/'inputs/sequence-accounting.tsv')}
assignments={r['accession']:r for r in read(R/'inputs/accession-assignments.tsv')}
strata={};members=[]
for r in read(R/'database/assembly-overlap.tsv'):
 acc=r['accession']
 status=('assembly_accession_detected_in_available_build_report' if r['exact_assembly_accession_in_library']=='1' or r['paired_assembly_accessions_in_library'] else 'sequence_accession_or_NZ_wrapper_detected_only' if int(r['exact_contigs']) or int(r['nz_wrapper_contigs']) else 'no_detected_accession_match_independence_unknown')
 type_status='study_type_anchor' if assignments[acc]['formal_type_names'] else 'not_study_type_anchor'
 strata[acc]={'database_accession_evidence':status,'type_anchor_role':type_status}
 members.append(dict(accession=acc,**strata[acc],independence='not_established_for_any_stratum'))
assert set(strata)==set(account) and len(strata)==1070
write(R/'profile-details/overlap-stratum-membership.tsv',members)
counts=collections.defaultdict(collections.Counter);fractions=collections.defaultdict(lambda:collections.defaultdict(float))
for r in read(R/'summaries/taxid-profiles.tsv'):
 for kind,status in strata[r['accession']].items():
  for dimension,category in [('resolution',r['resolution']),('framework',r['framework_comparison'])]:
   key=r['input_mode'],r['confidence'],kind,status,dimension,category
   counts[key].update(sequences=int(r['sequences']),bases=int(r['bases']))
   fractions[key][r['accession']]+=float(r['represented_base_fraction'])
rows=[]
for key,c in sorted(counts.items()):
 mode,conf,kind,status,dimension,category=key
 accs=[acc for acc in account if strata[acc][kind]==status]
 base=sum(int(account[acc]['bases']) for acc in accs)
 seq=sum(int(account[acc]['contigs' if mode=='native' else 'fragments']) for acc in accs)
 rows.append(dict(input_mode=mode,confidence=conf,stratum_kind=kind,stratum=status,dimension=dimension,category=category,assembly_denominator=len(accs),sequence_denominator=seq,base_denominator=base,pooled_sequences=c['sequences'],pooled_represented_bases=c['bases'],pooled_sequence_fraction=c['sequences']/seq,pooled_base_fraction=c['bases']/base,mean_assembly_base_fraction=statistics.mean(fractions[key].get(acc,0) for acc in accs),independence='not_established_for_any_stratum'))
write(R/'profile-details/profile-summary-by-provenance-stratum.tsv',rows)
print(collections.Counter(r['database_accession_evidence'] for r in members))
print(collections.Counter(r['type_anchor_role'] for r in members))
