#!/usr/bin/env python3
"""Classify framework-relative label changes without inferring new genome assignments."""
import argparse,collections,csv,hashlib,json,re,time
from pathlib import Path
CATEGORIES=['unresolved_caution','biological_reassignment','framework_rank_change','nomenclature_only','added_resolution','retained_assignment']
INELIGIBLE='ineligible'

def read(path):
 with open(path) as handle:return list(csv.DictReader(handle,delimiter='\t'))
def write(path,rows,fields=None):
 with open(path,'w') as handle:
  w=csv.DictWriter(handle,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def parse_label(raw):
 """Preserve raw separately; parse explicit binomial/subspecies, not strain suffixes."""
 text=raw.strip().strip('"').strip()
 if not text:return {'status':'missing','species':'','subspecies':'','name':''}
 if re.search(r'\bsp\.|\bcomplex\b|\b(unclassified|unidentified|bacterium)\b',text,re.I):
  return {'status':'unspecified','species':'','subspecies':'','name':''}
 m=re.match(r'^([A-Z][a-z]+ [a-z][a-z-]+)(?: subsp\.? ([a-z][a-z-]+))?(?:\s|$)',text)
 if not m:return {'status':'unmapped_text','species':'','subspecies':'','name':''}
 species,sub=m.group(1),m.group(2) or ''
 return {'status':'named','species':species,'subspecies':sub,'name':species+(' subsp. '+sub if sub else '')}

def make_crosswalk(rows):
 species=collections.defaultdict(set);exact=collections.defaultdict(set);anchors=collections.defaultdict(set);direct=collections.defaultdict(set)
 for row in rows:
  target=row['project_species']
  if not target or target.endswith(' sp.'):continue
  for raw in row['formal_type_names'].split(';'):
   if not raw:continue
   parsed=parse_label(raw)
   if parsed['status']!='named':raise ValueError('Unparsed type name '+raw)
   species[parsed['species']].add(target)
   exact[parsed['name']].add((target,row['project_subspecies']))
   anchors[parsed['name']].add(row['accession'])
   if not parsed['subspecies'] or parsed['subspecies']==parsed['species'].split()[-1]:direct[parsed['species']].add(target)
 # Explicit project-only name: approved adelaidei core lacks a retained nomenclatural type.
 species['Enterobacter adelaidei'].add('Enterobacter adelaidei')
 exact['Enterobacter adelaidei'].add(('Enterobacter adelaidei',''))
 # Lelliottia named context is an explicit reviewed treatment; identity retention remains qualified.
 species['Lelliottia nimipressuralis'].add('Lelliottia nimipressuralis')
 exact['Lelliottia nimipressuralis'].add(('Lelliottia nimipressuralis',''))
 return species,exact,anchors,direct

def classify(row,raw,crosswalk,aliases=None):
 species,exact,anchors,direct=crosswalk;aliases=aliases or {}
 parsed=parse_label(raw);src=parsed['species'];target=row['project_species'];sub=row['project_subspecies'];flags=[]
 eligible=row['in_frozen_study']=='1' and row['assignment_status']!='excluded_or_restricted_analysis_scope'
 target_named=bool(target and not target.endswith(' sp.'))
 old_targets=species.get(src,set());substatus='not_explicitly_asserted';old_groups=exact.get(parsed['name'],set())
 # Aliases here require explicit same-circumscription/rank evidence; none inferred from GTDB names.
 if parsed['status']=='named' and src in aliases:
  old_targets=species.get(aliases[src],set())
 if not eligible:
  cat=INELIGIBLE;reason='outside_frozen_study' if row['in_frozen_study']!='1' else 'retained_QC_analysis_exclusion'
 elif not target_named:
  cat='unresolved_caution';reason='approved_unnamed_or_genus_qualified_outcome' if row['assignment_status']=='approved_review_treatment' else 'retained_context_without_reviewed_named_assignment'
 elif parsed['status'] in ('missing','unspecified'):
  cat='added_resolution';reason='source_has_no_explicit_species_assignment'
 elif not old_targets:
  cat='unresolved_caution';reason='old_named_label_has_no_validated_framework_crosswalk'
 elif target not in old_targets:
  cat='biological_reassignment';reason='approved_target_outside_old_species_groups_in_same_framework'
 elif src in aliases and aliases[src]==target:
  cat='nomenclature_only';reason='documented_same_rank_same_circumscription_alias'
 elif src!=target:
  cat='framework_rank_change';reason='old_type_group_transferred_or_merged_into_project_species'
 else:
  cat='retained_assignment';reason='same_species_label_supported_under_reviewed_framework'
 species_cat=cat;species_reason=reason
 if len(old_targets)>1:flags.append('old_species_spans_multiple_project_species')
 if target in old_targets and src!=target:flags.append('framework_name_or_rank_change')
 # Identify a broadened final species without pretending an unchanged binomial changed text.
 incoming={k for k,v in species.items() if target in v}
 if len(incoming)>1 and target_named:flags.append('target_has_revised_circumscription')
 if parsed['subspecies'] and eligible and target_named:
  if not old_groups:
   substatus='unresolved_old_subspecies_mapping';cat='unresolved_caution';reason='explicit_old_subspecies_lacks_type_crosswalk'
  elif target not in {x[0] for x in old_groups}:
   substatus='biological_reassignment_across_project_species';cat='biological_reassignment';reason='explicit_old_subspecies_type_group_differs_from_approved_species'
  elif sub and (target,sub) not in old_groups and all(x[1] for x in old_groups if x[0]==target):
   substatus='biological_reassignment_within_project_species';cat='biological_reassignment';reason='explicit_old_subspecies_type_group_differs_from_approved_subspecies'
  elif not sub:
   substatus='project_species_only_old_subspecies_not_reassessed';flags.append('old_subspecies_not_claimed_wrong_from_missing_project_subspecies')
  elif parsed['name']!=target+' subsp. '+sub:
   substatus='framework_rank_change'
   if cat=='retained_assignment':cat='framework_rank_change';reason='explicit_subspecies_transferred_or_merged_within_retained_species'
  else:substatus='retained_assignment'
 elif eligible and target_named and sub:flags.append('project_subspecies_resolution_not_asserted_by_source')
 if parsed['subspecies']:flags.append('source_explicit_subspecies')
 if not target_named:flags.append('no_named_species_resolution_not_an_open_author_vote')
 if row['formal_type_names']:flags.append('type_anchor_self_not_performance_validation')
 if row['nomenclature_and_boundary_flags']:flags.append('approved_qualifications_preserved')
 if row['assignment_status']=='source_name_retained_not_individually_adjudicated':flags.append('retention_not_individual_author_adjudication')
 if cat=='biological_reassignment':
  if row['assignment_status']!='approved_review_treatment':
   raise ValueError('Biological reassignment lacks explicit approval: '+row['accession'])
  flags.append('positive_approved_assignment_not_old_reference_failure_alone')
 strict=species_cat
 if eligible and parsed['status']=='named' and target_named and direct.get(src) and target not in direct[src]:strict='biological_reassignment'
 return dict(source_status=parsed['status'],source_species=src,source_subspecies=parsed['subspecies'],source_canonical_name=parsed['name'],old_framework_species=';'.join(sorted(old_targets)),old_exact_framework_groups=';'.join(a+(' subsp. '+b if b else '') for a,b in sorted(old_groups)),old_label_type_accessions=';'.join(sorted(anchors.get(parsed['name'],set()))),eligible=int(eligible),primary_category=cat,primary_reason=reason,species_only_category=species_cat,species_only_reason=species_reason,subspecies_outcome=substatus,affected_rank=('subspecies_within_project_species' if substatus=='biological_reassignment_within_project_species' else 'species_from_explicit_subspecies' if substatus=='biological_reassignment_across_project_species' else 'species'),strict_nominate_species_sensitivity=strict,reason_flags=';'.join(flags),species_display_changes=int(bool(target_named and src and src!=target)),named_species_added=int(target_named and parsed['status'] in ('missing','unspecified')))

def main():
 start=time.monotonic();p=argparse.ArgumentParser(description=__doc__);p.add_argument('--assignments',required=True);p.add_argument('--historical-release-cases',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
 rows=read(a.assignments);assert len(rows)==len({x['accession'] for x in rows});cw=make_crosswalk(rows);cross=[]
 for name,groups in sorted(cw[0].items()):
  cross.append(dict(source_species=name,project_species_groups=';'.join(sorted(groups)),crosswalk_basis='documented_parent_species_including_retained_subspecies_types' if len(groups)>1 else 'retained_type_to_approved_project_group' if any(k==name or k.startswith(name+' subsp. ') for k in cw[2]) else 'explicit_project_name_no_formal_type_claim',direct_nominate_groups=';'.join(sorted(cw[3].get(name,set())))))
 write(out/'species-framework-crosswalk.tsv',cross)
 exactrows=[dict(source_name=k,project_groups=';'.join(s+(' subsp. '+ss if ss else '') for s,ss in sorted(v)),type_accessions=';'.join(sorted(cw[2].get(k,set())))) for k,v in sorted(cw[1].items())]
 write(out/'exact-name-crosswalk.tsv',exactrows)
 records=[]
 for row in rows:
  for snapshot,column,date in [('historical','historical_ncbi','student_snapshot_retrieval_date_not_established'),('current_2026-09-23','current_ncbi',row['current_ncbi_date'])]:
   record=dict(row);record.update(snapshot=snapshot,original_label=row[column],source_snapshot_date=date,source_taxid='not_present_in_assignment_input');record.update(classify(row,row[column],cw));records.append(record)
 write(out/'all-records-categorised.tsv',records)
 study=[x for x in records if x['in_frozen_study']=='1'];assert len(study)==2148
 write(out/'study-categories.tsv',study)
 write(out/'eligible-accessions.tsv',[dict(accession=x['accession'],is_study_type=int(bool(x['formal_type_names'])),role=x['record_role']) for x in rows if x['in_frozen_study']=='1' and x['assignment_status']!='excluded_or_restricted_analysis_scope'])
 write(out/'excluded-accessions.tsv',[dict(accession=x['accession'],reason=x['historical_exclusion_context_status'],historical_label=x['historical_ncbi']) for x in rows if x['in_frozen_study']=='1' and x['assignment_status']=='excluded_or_restricted_analysis_scope'])
 write(out/'nonstudy-roles.tsv',[dict(accession=x['accession'],role=x['record_role'],status=x['assignment_status']) for x in rows if x['in_frozen_study']!='1'])
 counts=[];members=[]
 for snap in ['historical','current_2026-09-23']:
  s=[x for x in study if x['snapshot']==snap]
  cohorts={'frozen_study_1074':s,'eligible_study_1070':[x for x in s if x['eligible']], 'eligible_non_type':[x for x in s if x['eligible'] and not x['formal_type_names']], 'assessable_previously_named':[x for x in s if x['eligible'] and x['source_status']=='named' and x['primary_category']!='unresolved_caution']}
  for cohort,subset in cohorts.items():
   for x in subset:members.append(dict(snapshot=snap,denominator=cohort,accession=x['accession']))
   for view in ['primary_category','species_only_category','strict_nominate_species_sensitivity']:
    ctr=collections.Counter(x[view] for x in subset)
    for category in [INELIGIBLE]+CATEGORIES:
     n=ctr[category];counts.append(dict(snapshot=snap,denominator_name=cohort,view=view,category=category,numerator=n,denominator=len(subset),percent=100*n/len(subset) if subset else '',rule='protocol_v1_2026-10-03'))
 write(out/'category-counts.tsv',counts);write(out/'denominator-membership.tsv',members)
 review=[x for x in study if x['primary_category']!=x['species_only_category'] or x['strict_nominate_species_sensitivity']!=x['species_only_category'] or x['primary_category']=='unresolved_caution']
 write(out/'implementation-review-cases.tsv',review,list(study[0]))
 transitions=[]
 for snap in ['historical','current_2026-09-23']:
  c=collections.Counter((x['source_species'] or '['+x['source_status']+']',x['project_species'] or '[no named assignment]',x['primary_category'],x['affected_rank']) for x in study if x['snapshot']==snap)
  for (old,new,category,rank),n in sorted(c.items()):transitions.append(dict(snapshot=snap,old_species=old,new_species=new,category=category,affected_rank=rank,count=n))
 write(out/'label-transitions.tsv',transitions)
 # Compact evidence extracts preserve approval, rather than using the best ANI hit as a classifier.
 fields=['accession','original_label','project_species','project_subspecies','primary_category','primary_reason','affected_rank','subspecies_outcome','old_framework_species','old_exact_framework_groups','old_label_type_accessions','assignment_status','decision_authority','rule_id','nomenclature_and_boundary_flags']
 historical=[x for x in study if x['snapshot']=='historical' and x['eligible']]
 for filename,subset in [('biological-reassignments.tsv',[x for x in historical if x['primary_category']=='biological_reassignment']),('subspecies-reassignments.tsv',[x for x in historical if x['primary_category']=='biological_reassignment' and x['source_subspecies']]),('caution-cases.tsv',[x for x in historical if x['primary_category']=='unresolved_caution'])]:
  write(out/filename,[{k:x[k] for k in fields} for x in subset],fields)
 flagcounts=collections.Counter(flag for x in historical for flag in x['reason_flags'].split(';') if flag)
 write(out/'reason-flag-counts.tsv',[dict(flag=k,count=v,denominator=len(historical)) for k,v in sorted(flagcounts.items())])
 # Practical same-release comparison uses source strings unchanged. Join only the saved audit fields.
 release={x['accession']:x for x in read(a.historical_release_cases)};practical=[]
 for x in rows:
  if x['in_frozen_study']!='1':continue
  saved=release[x['accession']];assert saved['historical_kraken_chromosome']==x['historical_kraken_chromosome'] and saved['historical_kraken_whole']==x['historical_kraken_whole']
  ch=x['historical_kraken_chromosome'];wh=x['historical_kraken_whole'];missing=not ch or not wh or any(v.lower() in ('unclassified','unassigned','unknown') for v in [ch,wh])
  practical.append(dict(accession=x['accession'],eligible=int(x['assignment_status']!='excluded_or_restricted_analysis_scope'),chromosome_call=ch,whole_call=wh,comparison_status='missing_or_unclassified' if missing else 'disagreement' if ch!=wh else 'agreement',release='GTDB_R220_historical_database',r220_accession_match=saved['r220_accession_match'],r220_record_accession=saved['r220_record_accession'],r220_accession_gtdb_species=saved['r220_accession_gtdb_species'],chromosome_vs_r220=saved['chromosome_vs_r220_accession'],whole_vs_r220=saved['whole_vs_r220_accession']))
 write(out/'practical-disagreement.tsv',practical)
 practical_counts=[]
 for scope,subset in [('frozen_study_1074',practical),('eligible_study_1070',[x for x in practical if x['eligible']])]:
  for status in ['agreement','disagreement','missing_or_unclassified']:
   n=sum(x['comparison_status']==status for x in subset);practical_counts.append(dict(denominator_name=scope,status=status,numerator=n,denominator=len(subset),percent=100*n/len(subset)))
 write(out/'practical-counts.tsv',practical_counts)
 inputs=[a.assignments,a.historical_release_cases,__file__];write(out/'input-hashes.tsv',[dict(path=str(p),sha256=sha(p)) for p in inputs])
 summary=dict(study_records=len(rows)-sum(x['in_frozen_study']!='1' for x in rows),all_audit_records=len(rows),category_rows=len(records),elapsed_seconds=round(time.monotonic()-start,3),historical_eligible_counts=dict(collections.Counter(x['primary_category'] for x in study if x['snapshot']=='historical' and x['eligible'])),current_eligible_counts=dict(collections.Counter(x['primary_category'] for x in study if x['snapshot']=='current_2026-09-23' and x['eligible'])),practical_counts=practical_counts)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
