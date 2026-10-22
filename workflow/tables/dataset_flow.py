"""Count/selection routine from the frozen manuscript renderer; no plotting or private receipts."""
def dataset_flow():
    a = rows(PATHS['accessions'])
    study = {x['accession'] for x in a if x['in_frozen_study'] == '1'}
    roles = Counter((x['record_role'] for x in a if x['in_frozen_study'] == '0'))
    restricted = {x['accession'] for x in a if x['in_frozen_study'] == '1' and x['assignment_status'] == 'excluded_or_restricted_analysis_scope'}
    ref = rows(PATHS['reference_panel'])
    refset = {x['accession'] for x in ref}
    split = rows(PATHS['development_split'])
    native = rows(PATHS['native_execution'])
    studyquery = {x['genome_id'] for x in native if x['layer'] == 'layer1'}
    publicquery = {x['genome_id'] for x in native if x['layer'] == 'layer3'}
    assert studyquery == study - restricted - refset
    assert len(refset & study) == 124 and len(refset - study) == 16
    assert len(split) == len(ref) + len(studyquery) == 1086
    members = rows(PATHS['genomovar_members'])
    memberids = {x['accession'] for x in members}
    assert len(memberids) == 173 and len(memberids & study) == 172 and (len(memberids - study) == 1)
    counts = {'genomovar_members': len(memberids), 'genomovar_study_members': len(memberids & study), 'genomovar_additional_types': len(memberids - study), 'initial_collection': 1079, 'initial_exclusions': roles['task01_public_exclusion'], 'study': len(study), 'final_quality_exclusions': len(restricted), 'classification': len(study - restricted), 'catalogue': len(a), 'additional_types': roles['type_anchor_nonstudy'], 'outgroups': roles['outgroup'], 'diagnostic_addition': roles['diagnostic_only_addition'], 'reference_selection': len(split), 'primary_references': len(ref), 'study_references': len(refset & study), 'study_queries': len(studyquery), 'development': sum((x['partition'] == 'development' for x in split)), 'internal_reserve': sum((x['partition'] != 'development' for x in split)), 'public_pool': len(rows(PATHS['public_pool'])), 'confirmation': len(rows(PATHS['confirmation'])), 'linked_hold': len(rows(PATHS['linked_hold'])), 'public_evaluation': len(publicquery), 'public_source_comparison': len(rows(PATHS['public_truth'])), 'public_exposed': len(rows(PATHS['public_held_truth'])), 'native_total': len(native), 'expert_subset': json.loads(PATHS['native_preparation'].read_text())['expert']}
    assert counts['initial_collection'] == counts['study'] + counts['initial_exclusions']
    assert counts['catalogue'] == counts['study'] + sum(roles.values())
    assert counts['public_pool'] == counts['confirmation'] + counts['linked_hold'] + counts['public_evaluation'] == 6424
    assert counts['public_evaluation'] == len(rows(PATHS['public_reserve'])) == counts['public_source_comparison'] + counts['public_exposed'] == 4857
    assert counts['native_total'] == counts['study_queries'] + counts['public_evaluation'] == 5803
    transition = rows(PATHS['initial_transitions'])
    assert int(next((x['count'] for x in transition if x['rule'] == 'public_discovery:entities'))) == 1079
    with (tab / 'S01-dataset-flow.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['dataset', 'assemblies_or_records'])
        w.writerows(counts.items())
    fields = ['accession', 'project_species', 'project_subspecies', 'type_anchor', 'source_role', 'selection_reason', 'completeness', 'contamination', 'quality_flags', 'whole_sha256']
    with (tab / 'S08-primary-references.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(ref)
