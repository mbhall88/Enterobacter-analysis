#!/usr/bin/env python3
"""Final independent counts, conservation, commands, selected-pair and resource checks."""
import csv,json,collections
from datetime import datetime
from pathlib import Path
import task12_scope_s1_s2 as s
R=s.RUN
exc=s.rows(s.AUD/'excluded-records.tsv');r=s.rows(R/'s1-records.tsv');c=s.rows(R/'s2-contigs.tsv');a=s.rows(R/'s2-anchors.tsv');checks={}
checks['excluded_unique_complete']=len(r)==len(exc)==len({(x['accession'],x['record']) for x in r})==2839
checks['excluded_bases_conserved']=sum(int(x['bp']) for x in r)==sum(int(x['bp']) for x in exc)==251178419
for x in r:
 if x['status']=='screened':assert 0<=int(x['shared_chromosome_plasmid_bp'])<=min(int(x['chromosome_aligned_bp']),int(x['plasmid_aligned_bp']))<=int(x['bp'])
checks['label_overlap_bounds']=True
checks['draft_anchors_complete']=len(a)==16 and len(c)==891 and len({(x['accession'],x['record']) for x in c})==891
checks['draft_bases_conserved']=sum(int(x['bp']) for x in c)==sum(int(x['draft_bp']) for x in a)
logs=list((R/'alignments').glob('*/commands.json'));checks['screen_commands_success']=len(logs)==762 and all(x['exit_status']==0 for p in logs for x in json.loads(p.read_text()))
checks['screen_commands_expected']=sum(len(json.loads(p.read_text())) for p in logs)==748*6+14*3
jobs=s.rows(R/'s2-anim-jobs.tsv');checks['sensitivity_jobs_complete']=len(jobs)==18 and all(json.loads((Path(j['reuse_source'] if j['reused']=='True' else j['output'])/'command.json').read_text())['exit_code']==0 for j in jobs)
checks['sensitivity_directions_complete']=len(s.rows(R/'s2-anim-raw-directions.tsv'))==36
checks['original_hashes_still_match']=all(s.hashfile(x['original_path'])==x['original_sha256'] for x in s.rows(R/'s2-derived-inputs.tsv'))
tips={x['accession'] for x in s.rows(s.AUD/'tip-scope-accounting.tsv')};census=s.rows(R/'s2-pair-selection-census.tsv');pairs={(x['type_accession'],x['query_accession']) for x in census};missing=[dict(type_accession=x['accession'],query_accession=q,status='frozen_tip_absent_from_existing_skani_ledger') for x in a for q in sorted(tips-{x['accession']}) if (x['accession'],q) not in pairs];s.write(R/'s2-unmeasured-frozen-pairs.tsv',missing,['type_accession','query_accession','status'])
checks['pair_census_unique']=len(pairs)==len(census)==17424
checks['frozen_pair_census_complete']=sum(x['query_accession'] in tips for x in census)+len(missing)==16*1088
flag=[x for x in r if x['classification']=='chromosome_like'];summary=dict(s1_records=len(r),s1_bp=sum(int(x['bp']) for x in r),s1_status=collections.Counter(x['status'] for x in r),s1_class=collections.Counter(x['classification'] for x in r),s1_flagged_assemblies=len({x['accession'] for x in flag}),s1_flagged_bp=sum(int(x['bp']) for x in flag),s1_co_resident_90pct_flags=[x['accession'] for x in flag if int(x['self_chromosome_aligned_bp'])/int(x['bp'])>=.9],s2_anchors=len(a),s2_contigs=len(c),s2_bp=sum(int(x['bp']) for x in c),s2_class=collections.Counter(x['classification'] for x in c),s2_pair_dispositions=collections.Counter(x['disposition'] for x in census),s2_missing_frozen_pairs=len(missing),checks=checks)
# Scheduler allocation peak from job rows, not requested CPU flags.
with open(R/'scheduler/anim-sacct.tsv') as f:acc=[x for x in csv.DictReader(f,delimiter='|') if '.' not in x['JobID']]
checks['scheduler_anim_complete']=len(acc)==17 and all(x['State']=='COMPLETED' and x['ExitCode']=='0:0' for x in acc)
# Slurm configured display dates use HH:MM:SS DD/MM/YY locally.
parse=lambda v:datetime.strptime(v,'%H:%M:%S %d/%m/%y') if ' ' in v else datetime.fromisoformat(v)
events=[]
for x in acc:events.extend([(parse(x['Start']),int(x['AllocCPUS']),1),(parse(x['End']),-int(x['AllocCPUS']),-1)])
cpu=n=peak=peakjobs=0
for dt,delta,jobdelta in sorted(events,key=lambda z:(z[0],z[1])):
 cpu+=delta;n+=jobdelta;peak=max(peak,cpu);peakjobs=max(peakjobs,n)
summary.update(anim_peak_actual_cpus=peak,anim_peak_jobs=peakjobs,anim_job_elapsed_seconds=[sum(int(n)*m for n,m in zip(x['Elapsed'].split(':'),[3600,60,1])) for x in acc])
assert all(checks.values()),checks
(R/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
