"""Reconcile saved cohort manifests and supply Figure 1 source data; no sequences opened."""
from pathlib import Path
from collections import Counter
import csv,json
from prepare import ROOT,RUN,sha,rows

DESIGN=ROOT/'results/tool-planning/2026-10-04/balanced-marker-discovery-astra'
PUBLIC=ROOT/'results/tool-planning/2026-10-06/evaluation-freeze-design-astra'
EVAL=ROOT/'results/tool-planning/2026-10-07/validation-runs-astra'
PATHS={
 'accessions':ROOT/'results/ani-cases/2026-10-03/task12-astra-cp102-approved-final-reconciliation/accession-assignments.tsv',
 'initial_transitions':ROOT/'results/cohort/2026-09-16/task01-final-manifests-opus/count-transitions.tsv',
 'reference_panel':DESIGN/'design-v2/discovery-panel.tsv',
 'development_split':DESIGN/'design-v2/split-manifest.tsv',
 'reference_methods':DESIGN/'manuscript-contribution.md',
 'reference_selection_code':DESIGN/'make_benchmark.snapshot.py',
 'genomovar_members':ROOT/'results/ani-cases/2026-10-02/task12-astra-genomovars-reporting/outputs/memberships.tsv',
 'public_pool':ROOT/'results/tool-planning/2026-10-06/conservative-boundary-validation-astra/later-pool-lossless/later-validation-metadata-pool-lossless.tsv',
 'confirmation':PUBLIC/'reserve-protected-alternative/outputs/confirmation-manifest.tsv',
 'linked_hold':PUBLIC/'reserve-protected-alternative/outputs/confirmation-linked-hold.tsv',
 'public_reserve':PUBLIC/'reserve-protected-alternative/outputs/remaining-reserve-metadata.tsv',
 'public_methods':PUBLIC/'protocol.md',
 'native_execution':EVAL/'inputs/execution-manifest.tsv',
 'study_truth':EVAL/'inputs/layer1-truth.tsv',
 'public_truth':EVAL/'inputs/layer3-truth.tsv',
 'public_held_truth':EVAL/'inputs/layer3-held-truth.tsv',
 'native_preparation':EVAL/'inputs/preparation.json',
 'evaluation_methods':ROOT/'results/tool-planning/2026-10-08/validation-adjudications-astra/validation-report.md',
 'ani_commands':ROOT/'results/ani-screen/2026-09-17/task04-opus-type-anchored/commands.sh',
 'ani_methods':ROOT/'results/ani-screen/2026-09-17/task04-opus-type-anchored/README.md',
}

def build():
    a=rows(PATHS['accessions']); study={x['accession'] for x in a if x['in_frozen_study']=='1'}
    roles=Counter(x['record_role'] for x in a if x['in_frozen_study']=='0')
    restricted={x['accession'] for x in a if x['in_frozen_study']=='1' and x['assignment_status']=='excluded_or_restricted_analysis_scope'}
    ref=rows(PATHS['reference_panel']); refset={x['accession'] for x in ref}; split=rows(PATHS['development_split'])
    native=rows(PATHS['native_execution']); studyquery={x['genome_id'] for x in native if x['layer']=='layer1'}
    publicquery={x['genome_id'] for x in native if x['layer']=='layer3'}
    assert studyquery==(study-restricted-refset)
    assert len(refset & study)==124 and len(refset-study)==16
    assert len(split)==len(ref)+len(studyquery)==1086
    members=rows(PATHS['genomovar_members']); memberids={x['accession'] for x in members}
    assert len(memberids)==173 and len(memberids & study)==172 and len(memberids-study)==1
    counts={'genomovar_members':len(memberids),'genomovar_study_members':len(memberids & study),'genomovar_additional_types':len(memberids-study),'initial_collection':1079,'initial_exclusions':roles['task01_public_exclusion'],
      'study':len(study),'final_quality_exclusions':len(restricted),'classification':len(study-restricted),
      'catalogue':len(a),'additional_types':roles['type_anchor_nonstudy'],'outgroups':roles['outgroup'],
      'diagnostic_addition':roles['diagnostic_only_addition'],'reference_selection':len(split),
      'primary_references':len(ref),'study_references':len(refset & study),'study_queries':len(studyquery),
      'development':sum(x['partition']=='development' for x in split),
      'internal_reserve':sum(x['partition']!='development' for x in split),
      'public_pool':len(rows(PATHS['public_pool'])),'confirmation':len(rows(PATHS['confirmation'])),
      'linked_hold':len(rows(PATHS['linked_hold'])),'public_evaluation':len(publicquery),
      'public_source_comparison':len(rows(PATHS['public_truth'])),'public_exposed':len(rows(PATHS['public_held_truth'])),
      'native_total':len(native),'expert_subset':json.loads(PATHS['native_preparation'].read_text())['expert']}
    assert counts['initial_collection']==counts['study']+counts['initial_exclusions']
    assert counts['catalogue']==counts['study']+sum(roles.values())
    assert counts['public_pool']==counts['confirmation']+counts['linked_hold']+counts['public_evaluation']==6424
    assert counts['public_evaluation']==len(rows(PATHS['public_reserve']))==counts['public_source_comparison']+counts['public_exposed']==4857
    assert counts['native_total']==counts['study_queries']+counts['public_evaluation']==5803
    transition=rows(PATHS['initial_transitions'])
    assert int(next(x['count'] for x in transition if x['rule']=='public_discovery:entities'))==1079
    out=RUN/'outputs';tab=out/'tables'
    receipt={'counts':counts,'sources':{k:{'path':str(p),'sha256':sha(p)} for k,p in PATHS.items() if p.is_file()},
             'interpretation':'Saved accession/role accounting; no independent-cohort sum across overlapping analyses',
             'quality_exclusions':sorted(restricted),'study_reference_intersection':sorted(study & refset)}
    (out/'dataset-flow.json').write_text(json.dumps(receipt,indent=2)+'\n')
    with (tab/'S01-dataset-flow.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['dataset','assemblies_or_records'])
        w.writerows(counts.items())
    fields=['accession','project_species','project_subspecies','type_anchor','source_role','selection_reason','completeness','contamination','quality_flags','whole_sha256']
    with (tab/'S08-primary-references.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(ref)
    print(json.dumps(counts,indent=2))

if __name__=='__main__':build()
