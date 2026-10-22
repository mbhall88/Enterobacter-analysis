#!/usr/bin/env python3
"""Independent accounting checks; does not import the categorisation implementation."""
import collections,csv,hashlib,json
from pathlib import Path
run=Path(__file__).resolve().parent
root=run.parents[3]
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
checks=[]
def check(condition, name):
 assert condition,name
 checks.append(name)
sources=read(root/'results/ani-cases/2026-10-03/task12-astra-cp102-approved-final-reconciliation/accession-assignments.tsv')
original={r['accession']:r for r in sources}
rows=read(run/'all-records-categorised.tsv')
check(len(rows)==2196 and len(original)==1098,'one row per accession per source snapshot')
check(len({(r['snapshot'],r['accession']) for r in rows})==2196,'unique snapshot-accession keys')
check(all(all(r[k]==v for k,v in original[r['accession']].items()) for r in rows),'every original field including decisions, boundary flags, labels, paths and hashes preserved')
check(all(r['original_label']==r['historical_ncbi' if r['snapshot']=='historical' else 'current_ncbi'] for r in rows),'dated source labels not mixed')
expected_excluded={'GCA_002356475.1','GCA_006385655.1','GCA_041890495.1','GCA_050242595.1'}
study={r['accession'] for r in sources if r['in_frozen_study']=='1'}
check(len(study)==1074,'frozen study membership')
check({r['accession'] for r in read(run/'excluded-accessions.tsv')}==expected_excluded,'four later exclusions separately accounted')
check({r['accession'] for r in read(run/'eligible-accessions.tsv')}==study-expected_excluded,'exact eligible set without double subtraction')
check(len(read(run/'nonstudy-roles.tsv'))==24,'24 reference, outgroup, diagnostic or earlier exclusion records outside frozen denominator')
check(all(r['accession'].startswith('GCA_') for r in sources),'public assembly accessions only')
h=[r for r in rows if r['snapshot']=='historical' and r['eligible']=='1']
c=[r for r in rows if r['snapshot']=='current_2026-09-23' and r['eligible']=='1']
check(len(h)==len(c)==1070,'both eligible denominators')
check(sum(bool(r['formal_type_names']) for r in h)==17,'17 included study type anchors explicitly flagged')
check(all(r['assignment_status']=='approved_review_treatment' for r in h if r['primary_category']=='biological_reassignment'),'every reassignment has an approved placement')
check(all(r['primary_category']=='unresolved_caution' for r in h if not r['project_species'] or r['project_species'].endswith(' sp.')),'all final unnamed/context cases retained as caution')
check(sum(r['historical_ncbi']!=r['current_ncbi'] for r in sources if r['in_frozen_study']=='1')==1,'one source snapshot change; final categories stable independently checked')
check({r['accession']:r['primary_category'] for r in h}=={r['accession']:r['primary_category'] for r in c},'per-accession category concordance across snapshots')
# Independent tally via csvtk output, with no use of the implementation classifier.
frequency=read(run/'logs/csvtk-frequency.tsv')
observed=collections.Counter((r['snapshot'],r['eligible'],r['primary_category']) for r in rows if r['in_frozen_study']=='1')
check({(r['snapshot'],r['eligible'],r['primary_category']):int(r['frequency']) for r in frequency}==dict(observed),'csvtk independently reproduces all category frequencies')
counts=read(run/'category-counts.tsv'); members=read(run/'denominator-membership.tsv')
bykey={(r['snapshot'],r['accession']):r for r in rows}
for (snap,denom,view),rr in __import__('itertools').groupby(counts,key=lambda r:(r['snapshot'],r['denominator_name'],r['view'])):
 rr=list(rr); subset=[bykey[(snap,m['accession'])] for m in members if m['snapshot']==snap and m['denominator']==denom]
 independent=collections.Counter(r[view] for r in subset)
 check(len(subset)==len({r['accession'] for r in subset}) and all(int(r['denominator'])==len(subset) and int(r['numerator'])==independent[r['category']] and abs(float(r['percent'])-100*independent[r['category']]/len(subset))<1e-10 for r in rr) and sum(int(r['numerator']) for r in rr)==len(subset),f'count/percent/exhaustive partition: {snap}/{denom}/{view}')
# Explicit-subspecies evidence cannot be lost in binomial-only accounting.
subbio=[r for r in h if r['primary_category']!=r['species_only_category']]
check(len(subbio)==29 and all(r['source_subspecies'] and r['primary_category']=='biological_reassignment' and r['species_only_category']=='framework_rank_change' for r in subbio),'29 explicit-subspecies reassignment distinctions')
check(collections.Counter(r['subspecies_outcome'] for r in subbio)=={'biological_reassignment_within_project_species':28,'biological_reassignment_across_project_species':1},'28 within-species and one cross-species subgroup case')
check(sum(r['source_canonical_name']=='Enterobacter hormaechei' and r['project_species']=='Enterobacter xiangfangensis' for r in h)==434,'434 unqualified historical hormaechei treated conservatively as framework changes')
check(all(r['primary_category']=='framework_rank_change' for r in h if r['source_canonical_name']=='Enterobacter hormaechei' and r['project_species']=='Enterobacter xiangfangensis'),'no narrow-type-only inflation for bare hormaechei labels')
practical=read(run/'practical-disagreement.tsv')
check(len(practical)==1074 and {r['accession'] for r in practical}==study,'exact historical practical-comparison membership')
check(all(r['comparison_status']==('disagreement' if r['chromosome_call']!=r['whole_call'] else 'agreement') for r in practical),'independent raw same-release scope comparison')
check(sum(r['comparison_status']=='disagreement' for r in practical)==254 and sum(r['comparison_status']=='disagreement' and r['eligible']=='1' for r in practical)==254,'254 disagreements retained in both denominators')
for r in read(run/'input-hashes.tsv'):
 p=Path(r['path']);p=p if p.is_absolute() else root/p
 check(hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],'hash '+r['path'])
print(json.dumps(dict(status='PASS',checks=len(checks),details=checks),indent=2))
