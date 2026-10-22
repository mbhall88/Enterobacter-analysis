#!/usr/bin/env python3
"""Reconcile recorded approvals with the frozen public roster; do not infer approvals."""
import argparse,collections,csv,hashlib,json,re,time
from pathlib import Path
STUDY='manifests/study-candidates.tsv'
FULL='manifests/full.tsv'
CASES='results/ani-cases/2026-09-23/task12-astra-evidence/case-register/enriched-v3-complete/cases.tsv'
SEQUENCES='results/ani-cases/2026-09-23/task12-astra-evidence/inputs/resolved/genomes.tsv'
POPULATIONS='results/ani-cases/2026-10-02/task12-astra-genomovars-reporting/outputs/memberships.tsv'
DIAGNOSTICS='results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded/inputs/panel-primary.tsv'
INVENTORY='results/project-review/2026-09-29/task18-claude-pre-leah-opus55/checks/WS-I-inventory/accession-coverage.tsv'

def read(p):
 with open(p) as handle:return list(csv.DictReader(handle,delimiter='\t'))
def keyed(rows,key):
 d={}
 for r in rows:
  if r[key] in d:raise ValueError('Duplicate key '+r[key])
  d[r[key]]=r
 return d

def adjudications(rows):
 """Conflicting rules are an error, not a last-write-wins resolution."""
 out={}
 for r in rows:
  if r['accession'] in out:raise ValueError('Multiple assignment rules for '+r['accession'])
  out[r['accession']]=r
 return out

def exact_species(label):
 """Only a fully specified binomial; no 'complex sp.' or free-text inference."""
 m=re.fullmatch(r'(Enterobacter|Lelliottia) ([a-z][a-z]+)(?: subsp\. [a-z]+)?',label.strip())
 return ' '.join(m.groups()) if m else ''

def write(p,rows,fields=None):
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def reconcile(study,full,cases,seq,inventory,rules):
 outputs=[]
 universe=set(study)|set(inventory)
 unknown=set(rules)-universe
 if unknown:raise ValueError('Rules outside frozen audit universe: '+str(sorted(unknown)))
 for a in sorted(universe):
  s=study.get(a,{});f=full.get(a,{});c=cases.get(a,{});i=inventory.get(a,{});q=seq.get(a,{});r=rules.get(a,{})
  ins=a in study
  role='study' if ins else i.get('row_scope','reference_or_context')
  status='';species='';basis='';candidate=''
  if r:
   status='approved_review_treatment';species=r['project_species'];basis=r['decision_url']
   if 'context;' in r['qualification']: role='study_context' if ins else role
  elif i.get('coverage_class')=='excluded':
   status='excluded_or_restricted_analysis_scope';basis=i['exclusion_context_status']
  elif i.get('coverage_class')=='context_or_restricted':
   status='retained_context_restriction';basis=i['exclusion_context_status'] or f.get('decision_record','declared outgroup')
   if ins:role='study_context'
  elif i.get('is_type_anchor') in ('1','True','true'):
   status='retained_formal_type_anchor';basis='frozen type identity; no new naming inference'
   names={exact_species(n) for n in i.get('formal_type_names','').split(';')} - {''}
   if len(names)==1:species=names.pop()
  else:
   # A candidate is evidence, not an assignment. No highest-ANI auto-approval.
   pt=[t for t in c.get('passing_type_accessions','').split(';') if t]
   candidate_names=set()
   for t in pt:
    if t in rules and rules[t]['project_species']: candidate_names.add(rules[t]['project_species'])
    else:
     names=cases.get(t,{}).get('formal_type_names','').split(';')
     candidate_names.update(exact_species(n) for n in names if exact_species(n))
   candidate=';'.join(sorted(candidate_names))
   old=exact_species(c.get('current_ncbi',''))
   if len(candidate_names)==1 and old in candidate_names:
    status='source_name_retained_not_individually_adjudicated';species=old
    basis='unchanged current NCBI binomial; compatible passing-type evidence; no new author decision'
   else:
    status='assignment_review_required';basis='no explicit accession treatment located; candidate is not approved'
  rowsource=c.get('current_ncbi','')
  outputs.append(dict(accession=a,in_frozen_study=int(ins),record_role=role,assignment_status=status,
   project_species=species,project_subspecies=r.get('project_subspecies',''),review_lineage=r.get('lineage',''),
   nomenclature_and_boundary_flags=r.get('qualification',''),decision_authority=r.get('authority',''),
   decision_or_retention_basis=basis,rule_id=r.get('rule_id',''),rule_membership_source=r.get('evidence_table',''),
   rule_source_sha256=r.get('evidence_sha256',''),evidence_candidate_species=candidate,
   historical_ncbi=s.get('historical_ncbi',c.get('historical_ncbi','')),current_ncbi=rowsource,current_ncbi_date=c.get('current_ncbi_date',''),
   historical_kraken_chromosome=c.get('historical_kraken_chromosome',''),historical_kraken_whole=c.get('historical_kraken_whole',''),
   legacy_revised_species=s.get('legacy_revised_species',c.get('legacy_revised_species','')),legacy_revised_subspecies=s.get('legacy_revised_subspecies',c.get('legacy_revised_subspecies','')),
   gtdb_r232_species=c.get('gtdb_r232_species',''),formal_type_names=c.get('formal_type_names',i.get('formal_type_names','')),
   frozen_ingroup_status=f.get('ingroup_status',c.get('ingroup_status','')),usage_restrictions=f.get('usage_restrictions',c.get('restriction','')),
   historical_exclusion_context_status=i.get('exclusion_context_status',''),quality_flags=c.get('quality_flags',''),
   best_type_accession=c.get('best_type_accession',''),best_type_skani_ani=c.get('best_type_ani',''),
   best_type_af_query=c.get('best_type_af_query',''),best_type_af_reference=c.get('best_type_af_reference',''),
   passing_type_accessions=c.get('passing_type_accessions',''),evidence_status=c.get('evidence_status',''),
   bac120_anchor_context=c.get('bac120_fasttree_retained_anchor_clade_types',''),verticall_anchor_context=c.get('verticall_fastme_retained_anchor_clade_types',''),mashtree_anchor_context=c.get('mashtree_consensus_retained_anchor_clade_types',''),
   baseline_ani_fasta=q.get('source',''),baseline_ani_fasta_sha256=q.get('sha256',''),
   frozen_tree_sequence_scope=f.get('sequence_scope',''),frozen_tree_fasta=f.get('source_fasta_path',''),frozen_tree_fasta_sha256=f.get('source_fasta_sha256',''),
   source_evidence_table=CASES,prior_coverage_class=i.get('coverage_class','')))
 return outputs

def main():
 start=time.monotonic();p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rules',required=True);p.add_argument('--output',required=True);a=p.parse_args();d=Path(a.output);d.mkdir(parents=True,exist_ok=True)
 study=keyed(read(STUDY),'accession_version');full=keyed(read(FULL),'accession_version');cases=keyed(read(CASES),'accession');seq=keyed(read(SEQUENCES),'accession');inv=keyed(read(INVENTORY),'accession');rules=adjudications(read(a.rules))
 diagnostics=keyed(read(DIAGNOSTICS),'accession')
 assert 'GCA_030510045.1' not in study
 inv['GCA_030510045.1']={'accession':'GCA_030510045.1','row_scope':'diagnostic_only_addition'}
 seq['GCA_030510045.1']={k:diagnostics['GCA_030510045.1'][k] for k in ['source','sha256']}
 rows=reconcile(study,full,cases,seq,inv,rules)
 populations=keyed(read(POPULATIONS),'accession')
 for r in rows:
  pop=populations.get(r['accession'],{})
  for key in ['genomovar_id','genomovar_rule','genomovar_version','operational_population','ST','ST_status']:
   r[key]=pop.get(key,'')
  r['population_source']=POPULATIONS if pop else ''
 assert len(study)==1074
 assert {r['accession'] for r in rows if r['in_frozen_study']}==set(study)
 assert all(re.fullmatch(r'GC[AF]_\d+\.\d+',r['accession']) for r in rows)
 for r in rows:
  if r['in_frozen_study']:assert r['historical_ncbi']==study[r['accession']]['historical_ncbi']
  if r['assignment_status']=='assignment_review_required':assert not r['project_species']
 write(d/'accession-assignments.tsv',rows)
 write(d/'remaining-assignment-review.tsv',[r for r in rows if r['assignment_status']=='assignment_review_required'],list(rows[0]))
 write(d/'retained-source-names.tsv',[r for r in rows if r['assignment_status']=='source_name_retained_not_individually_adjudicated'],list(rows[0]))
 write(d/'scope-and-exclusions.tsv',[r for r in rows if 'context' in r['assignment_status'] or 'excluded' in r['assignment_status'] or r['record_role']=='study_context'],list(rows[0]))
 src=[STUDY,FULL,CASES,SEQUENCES,INVENTORY,DIAGNOSTICS,POPULATIONS,a.rules]
 write(d/'reconciliation-input-hashes.tsv',[dict(path=x,sha256=hashlib.sha256(Path(x).read_bytes()).hexdigest()) for x in src])
 summary={'total_records':len(rows),'frozen_study':len(study),'by_status_all':dict(collections.Counter(r['assignment_status'] for r in rows)),
 'by_status_study':dict(collections.Counter(r['assignment_status'] for r in rows if r['in_frozen_study'])),
 'approved_rules':len(rules),'unknown_rule_accessions':0,'unresolved_rule_conflicts':0,'missing_study_rows':0,
 'seconds':round(time.monotonic()-start,3),'note':'Accounted-for is not synonymous with approved. No final reporting denominator or percentages computed.'}
 (d/'audit-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
