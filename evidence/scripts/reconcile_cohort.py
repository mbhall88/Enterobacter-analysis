#!/usr/bin/env python3
"""Reconcile discovery and validation identities into explicit keep/hold/exclude manifests.

Joins the completed FASTA audit to the discovery roster crosswalk and the legacy
curated metadata. Public entities are keyed by accession.version; clinical validation
entities are keyed by an explicit sample alias that retains collection date and
assembly method. Duplicate or conflicting keys are reported, never silently collapsed.

Deliberate non-inferences: identical sequence fingerprints do not establish that two
records are the same isolate, and a shared sample number does not establish that two
files are the same sample. Both produce provenance questions and a `hold`, not a merge
and not an automatic choice of the larger or newer assembly.

Accepted adjudications and verified file-gap records enter through explicit optional inputs.
Without them the reconciler behaves exactly as before. With them, a decision is applied only if
it validates: known cohort action, recognised coordinator review status, exact accession.version
entities that exist in the cohort, an unambiguous canonical record, unique input keys, and file
records whose paths, checksums, sequence hashes, lengths and named accessions still verify on
disk. Anything else stops the run rather than producing a quiet manifest.
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime, timezone
from pathlib import Path

VALIDATION_FILENAME = re.compile(
    r'^(?P<run_prefix>.+?)_(?P<collection_date>\d{8})_Sample(?P<sample_number>\d+)_(?P<assembly_method>.+)\.(?:fasta|fna)$')
ACCESSION = re.compile(r'^(?P<db>GC[AF])_(?P<numeric>\d+)\.(?P<version>\d+)$')

PUBLIC_FILE_ROLES = {'public_whole': 'whole', 'study_chromosome_candidate': 'chromosome',
                     'classifier_reference_copy': 'classifier_copy'}
VALIDATION_FILE_ROLES = {'validation_candidate': 'candidate', 'validation_quarantine': 'quarantine',
                         'unexpected_query_in_study_dir': 'misfiled_copy'}
LEGACY_EXCLUDE_MARKERS = {'exclude?', 'exclude'}

QUESTION_ID = re.compile(r'^(?P<prefix>PQ-[A-Z]+-)(?P<number>\d+)$')

ADJUDICATION_FIELDS = ('decision_id', 'question_ids', 'entities', 'proposed_relationship',
                       'confidence', 'proposed_cohort_action', 'canonical_record_if_accepted',
                       'reason', 'coordinator_review')
RESOLVED_FILE_FIELDS = ('accession', 'question_ids', 'whole_path', 'whole_origin', 'whole_sha256',
                        'chromosome_path', 'chromosome_origin', 'chromosome_file_sha256',
                        'chromosome_sequence_sha256', 'chromosome_accession', 'chromosome_length',
                        'chromosome_header', 'ncbi_genome_notes', 'validation')

# Cohort actions this script knows how to apply. An accepted decision carrying any other action
# stops the run: silently ignoring an accepted decision is the failure this input exists to avoid.
PAIR_ACTIONS = frozenset({'one_inference_unit', 'retain_alternate_provenance'})
EXCLUDE_ACTIONS = frozenset({'exclude_from_phylogeny_and_classifier_reference', 'retain_in_audit'})
REPRESENTATION_ACTIONS = frozenset({'one_assembly_representation', 'retain_pair_crosswalk'})
KNOWN_ACTIONS = PAIR_ACTIONS | EXCLUDE_ACTIONS | REPRESENTATION_ACTIONS

# Coordinator review states. `accepted_by_coordinator*` is a decision; `no_scientific_ambiguity*`
# records a verified bookkeeping fact that needed no decision. Both are applied. A withheld state
# is recorded and changes nothing. Any other string is a validation failure.
REVIEW_ACCEPTED = 'accepted_by_coordinator'
REVIEW_NO_ADJUDICATION = 'no_scientific_ambiguity'
REVIEW_WITHHELD = ('proposed', 'pending', 'deferred', 'rejected', 'under_review')

ARCHIVE_ORIGINS = frozenset({'student_archive', 'student_archive_existing'})
RETRIEVED_ORIGIN_PREFIXES = ('ncbi_download', 'derived_')
UNRESOLVED_TAXONOMY_NOTE = 'unverified source organism'

# Question-status ranks. A question is never erased; the ledger keeps the most advanced
# annotation in `status` and every annotation, with its source, in `status_detail`.
STATUS_RANK = {'open': 0, 'deferred': 1, 'partially_accepted': 2, 'partially_resolved': 2,
               'accepted': 3, 'resolved': 4}

MANIFEST_NAMES = ('study-candidates.tsv', 'validation-candidates.tsv', 'exclusions.tsv')

CAVEATS = [
    'Sequence fingerprints are reused from the FASTA audit; they ignore headers, case, wrapping and '
    'contig order but not rotation or reverse complementation.',
    'Identical sequence is a provenance question, not evidence that two records are the same isolate.',
    'A shared sample number is not an isolate identity; aliases retain collection date and assembly method.',
    'No duplicate pair is resolved by preferring the larger or newer assembly; assembly history is required.',
    'Held entities stay out of downstream accuracy denominators until adjudicated.',
    'Legacy species strings are historical curation, not a validated naming policy.',
    'Resolved file-gap records are additional sequence representations with their own origins and '
    'checksums; they were not part of the 3,807-file FASTA audit and are reported separately.',
    'An adjudicated duplicate-submission pair is one inference unit, not a deleted accession: the '
    'alternate submission stays in provenance and exclusions with a link to the selected entity.',
]


class ReconcileError(ValueError):
    """A validation failure that must stop the run rather than produce a quiet manifest."""


def digest_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read_tsv_with_fields(path):
    """Read a TSV and return its rows with the header actually present in the file."""
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader), list(reader.fieldnames or [])


def require_fields(fields, required, source, path):
    missing = [field for field in required if field not in set(fields)]
    if missing:
        raise ReconcileError(f'{source} {path} is missing required columns: ' + ', '.join(missing))


def require_unique(rows, key, source, path):
    repeated = sorted(value for value, count in Counter(row[key] for row in rows).items() if count > 1)
    if repeated:
        raise ReconcileError(f'{source} {path} repeats {key}: ' + '; '.join(repeated))


def expand_question_ids(text):
    """Expand a `A;B` list and an inclusive `PQ-PUB-009..PQ-PUB-060` range into explicit ids."""
    ids = []
    for token in (part.strip() for part in (text or '').split(';')):
        if not token:
            continue
        if '..' in token:
            start, _, end = token.partition('..')
            first, last = QUESTION_ID.match(start.strip()), QUESTION_ID.match(end.strip())
            if not first or not last or first['prefix'] != last['prefix']:
                raise ReconcileError(f'Unreadable question-id range: {token}')
            if int(last['number']) < int(first['number']):
                raise ReconcileError(f'Descending question-id range: {token}')
            width = len(first['number'])
            ids.extend(f"{first['prefix']}{number:0{width}d}"
                       for number in range(int(first['number']), int(last['number']) + 1))
        elif QUESTION_ID.match(token):
            ids.append(token)
        else:
            raise ReconcileError(f'Unreadable question id: {token}')
    return ids


def read_fasta_records(path):
    """Read a FASTA the way the audit does: ascii, uppercase, wrapping and blank lines ignored."""
    records, name, chunks = [], None, []
    for line in Path(path).read_bytes().decode('ascii').splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                records.append((name, ''.join(chunks).upper()))
            name, chunks = line[1:].strip(), []
        elif name is None:
            raise ReconcileError(f'Sequence before FASTA header: {path}')
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, ''.join(chunks).upper()))
    if not records or any(not sequence for _, sequence in records):
        raise ReconcileError(f'Empty FASTA or empty sequence: {path}')
    return records


def sequence_multiset_sha256(records):
    """The FASTA audit's order- and header-insensitive fingerprint, for files it never saw.

    This is a different quantity from the plain sequence SHA-256 the file-gap resolver recorded;
    both are kept so a derived chromosome can be compared with either table.
    """
    parts = sorted((len(sequence), hashlib.sha256(sequence.encode('ascii')).hexdigest())
                   for _, sequence in records)
    return hashlib.sha256(json.dumps(parts, separators=(',', ':')).encode()).hexdigest()


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n',
                                extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, '') for field in fields})


def parse_accession(accession):
    """Split an accession.version; the numeric core pairs GCA and GCF representations."""
    match = ACCESSION.match(accession or '')
    if not match:
        return None
    return dict(accession_version=accession, accession_db=match['db'],
                assembly_numeric=match['numeric'], accession_version_number=int(match['version']))


def parse_validation_alias(filename):
    """Stable clinical alias: run prefix, collection date, unpadded sample number, assembly method.

    Sample number alone is not an identity: the archive holds one sample number collected
    on two different dates, and collapsing to the number would merge them.
    """
    match = VALIDATION_FILENAME.match(filename)
    if not match:
        return None
    parts = match.groupdict()
    parts['sample_number'] = int(parts['sample_number'])
    parts['zero_padded'] = bool(re.match(r'0\d', match['sample_number']))
    parts['alias'] = '|'.join([parts['run_prefix'], parts['collection_date'],
                               f"Sample{parts['sample_number']}", parts['assembly_method']])
    return parts


def index_unique(rows, key_of, source, must_agree=()):
    """Index by key, recording every repeated key instead of retaining the first row."""
    grouped = OrderedDict()
    for row in rows:
        grouped.setdefault(key_of(row), []).append(row)
    index, conflicts = OrderedDict(), []
    for key, group in grouped.items():
        index[key] = group[0]
        if len(group) == 1:
            continue
        disagreeing = sorted(field for field in must_agree
                             if len({row.get(field, '') for row in group}) > 1)
        conflicts.append(dict(
            source=source, key=key, rows=len(group),
            conflict_kind='duplicate_key_conflicting_values' if disagreeing else 'duplicate_key_identical_values',
            disagreeing_fields=';'.join(disagreeing),
            values=';'.join(sorted({str(row.get(field, '')) for row in group for field in disagreeing})),
            detail=';'.join(sorted(str(row.get('path', row)) for row in group))))
    return index, conflicts


def choose_canonical(files):
    """Deterministic canonical file for one alias: unpadded filename first, then shortest path.

    This is a naming convention, not an assembly-quality judgement, and it is only applied
    where every file for the alias carries the same sequence fingerprint.
    """
    return sorted(files, key=lambda row: (row['alias_parts']['zero_padded'],
                                          len(row['filename']), row['path']))[0]


def join_field(rows, field, separator=';'):
    return separator.join(str(row.get(field, '')) for row in rows)


class Questions:
    """Deterministic provenance-question register; ids are stable for stable inputs."""

    def __init__(self):
        self.rows = []
        self._counts = Counter()

    def add(self, cohort, rule, severity, entity_ids, file_count, evidence, evidence_source,
            resolution_required):
        code = {'public_discovery': 'PUB', 'clinical_validation': 'VAL'}[cohort]
        self._counts[code] += 1
        question_id = f'PQ-{code}-{self._counts[code]:03d}'
        self.rows.append(dict(question_id=question_id, cohort=cohort, rule=rule, severity=severity,
                              entities=len(entity_ids), entity_ids=';'.join(entity_ids),
                              file_count=file_count, evidence=evidence,
                              evidence_source=evidence_source,
                              resolution_required=resolution_required,
                              status='open', status_source='', status_detail=''))
        return question_id

    def annotate(self, question_ids, status, source, detail, entity_scope=None):
        """Record a status annotation without erasing the question or an earlier annotation.

        `entity_scope` is the set of entities the decision or evidence actually covers. A question
        whose entity list is only partly covered becomes `partially_*`, never fully settled.
        """
        wanted = set(question_ids)
        index = {row['question_id']: row for row in self.rows}
        unknown = sorted(wanted - set(index))
        if unknown:
            raise ReconcileError(f'{source} references unknown provenance questions: '
                                 + ', '.join(unknown))
        touched = []
        for question_id in sorted(wanted):
            row = index[question_id]
            listed = [entity for entity in row['entity_ids'].split(';') if entity]
            if entity_scope is None:
                applied, covered = status, len(listed)
            else:
                covered = len([entity for entity in listed if entity in entity_scope])
                if not covered:
                    continue
                applied = status if covered == len(listed) else f'partially_{status}'
            if applied not in STATUS_RANK:
                raise ReconcileError(f'Unknown question status: {applied}')
            note = f'{applied} [{source}] {detail} ({covered}/{len(listed)} entities)'
            row['status_detail'] = f"{row['status_detail']} || {note}".strip(' |') \
                if row['status_detail'] else note
            if STATUS_RANK[applied] > STATUS_RANK[row['status']]:
                row['status'], row['status_source'] = applied, source
            touched.append(question_id)
        return touched


def load_skani_evidence(paths):
    """Map an unordered pair of FASTA basenames to a reusable ANI/aligned-fraction string."""
    evidence = {}
    for path in paths:
        for row in read_tsv(path):
            if 'Ref_file' not in row or 'Query_file' not in row:
                raise ValueError(f'Not a skani dist table: {path}')
            pair = frozenset((Path(row['Ref_file']).name, Path(row['Query_file']).name))
            if len(pair) != 2:
                continue
            evidence.setdefault(pair, []).append(
                f"skani ANI={row.get('ANI', '')} AF_ref={row.get('Align_fraction_ref', '')} "
                f"AF_query={row.get('Align_fraction_query', '')} [{Path(path).name}]")
    return evidence


def pair_evidence(evidence, files):
    found = []
    for i, left in enumerate(files):
        for right in files[i + 1:]:
            found.extend(evidence.get(frozenset((left['filename'], right['filename'])), []))
    return '; '.join(found)


def build_public(audit_rows, crosswalk, legacy, questions):
    roster, roster_conflicts = index_unique(
        crosswalk, lambda row: row['accession_version'], 'crosswalk',
        must_agree=('historical_ncbi', 'legacy_revised_species', 'legacy_revised_subspecies'))
    legacy_index, legacy_conflicts = index_unique(
        legacy, lambda row: row['genome'], 'legacy_metadata',
        must_agree=('revised_species', 'revised_subspecies'))
    conflicts = roster_conflicts + legacy_conflicts

    files = defaultdict(lambda: defaultdict(list))
    for row in audit_rows:
        slot = PUBLIC_FILE_ROLES.get(row['role'])
        if slot:
            files[row['accession_version']][slot].append(row)
    for accession, slots in files.items():
        for slot, group in slots.items():
            if len(group) > 1:
                conflicts.append(dict(source=f'audit:{slot}', key=accession, rows=len(group),
                                      conflict_kind='duplicate_accession_within_file_role',
                                      disagreeing_fields='path', values='',
                                      detail=join_field(group, 'path')))

    entities = OrderedDict()
    for accession in sorted(set(roster) | set(files)):
        meta = parse_accession(accession) or dict(accession_db='', assembly_numeric='',
                                                  accession_version_number='')
        record = roster.get(accession, {})
        chromosome = (files[accession].get('chromosome') or [None])[0]
        whole = (files[accession].get('whole') or [None])[0]
        copy_file = (files[accession].get('classifier_copy') or [None])[0]
        legacy_row = legacy_index.get(Path(chromosome['filename']).stem) if chromosome else None
        entity = dict(
            entity_id=accession, cohort='public_discovery', accession_version=accession,
            accession_db=meta['accession_db'], assembly_numeric=meta['assembly_numeric'],
            in_discovery_roster=bool(record),
            historical_ncbi=record.get('historical_ncbi', ''),
            historical_kraken_chromosome=record.get('historical_kraken_chromosome', ''),
            historical_kraken_whole=record.get('historical_kraken_whole', ''),
            legacy_revised_species=record.get('legacy_revised_species', ''),
            legacy_revised_subspecies=record.get('legacy_revised_subspecies', ''),
            crosswalk_review_flag=record.get('initial_review_flag', ''),
            legacy_metadata_genome_id=legacy_row['genome'] if legacy_row else '',
            galah_95_cluster=legacy_row['galah_95_cluster'] if legacy_row else '',
            galah_98_cluster=legacy_row['galah_98_cluster'] if legacy_row else '',
            chromosome_path=chromosome['path'] if chromosome else '',
            chromosome_sha256=chromosome['sha256'] if chromosome else '',
            chromosome_sequence_sha256=chromosome['sequence_multiset_sha256'] if chromosome else '',
            chromosome_contigs=chromosome['contigs'] if chromosome else '',
            chromosome_total_bp=chromosome['total_bp'] if chromosome else '',
            chromosome_qc_flags=chromosome['qc_flags'] if chromosome else '',
            whole_path=whole['path'] if whole else '',
            whole_sha256=whole['sha256'] if whole else '',
            whole_sequence_sha256=whole['sequence_multiset_sha256'] if whole else '',
            whole_contigs=whole['contigs'] if whole else '',
            whole_total_bp=whole['total_bp'] if whole else '',
            whole_qc_flags=whole['qc_flags'] if whole else '',
            classifier_copy_path=copy_file['path'] if copy_file else '',
            classifier_copy_matches_chromosome=bool(
                copy_file and chromosome and copy_file['sha256'] == chromosome['sha256']),
            files=sum(len(group) for group in files[accession].values()),
            resolved_whole_path='', resolved_whole_origin='', resolved_whole_sha256='',
            resolved_chromosome_path='', resolved_chromosome_origin='',
            resolved_chromosome_file_sha256='', resolved_chromosome_sequence_sha256='',
            resolved_chromosome_multiset_sha256='', resolved_chromosome_accession='',
            resolved_chromosome_length='', resolved_chromosome_header='',
            resolved_ncbi_genome_notes='', resolved_file_checks='', resolved_file_source='',
            entity_role='study_candidate', usage_restrictions=[], canonical_public_record='',
            selected_inference_entity_id='', alternate_entity_ids=[], adjudication_ids=[],
            disposition='', disposition_reasons=[], notes=[], provenance_question_ids=[])
        if record and record.get('legacy_revised_species', '').strip().lower() in LEGACY_EXCLUDE_MARKERS:
            entity['legacy_exclude_marker'] = True
        entities[accession] = entity

    for entity in entities.values():
        if not entity['in_discovery_roster']:
            entity['notes'].append('file_present_without_roster_row')
        if not entity['chromosome_path'] and not entity['whole_path']:
            entity['disposition_reasons'].append('no_local_sequence_for_roster_accession')
        elif not entity['chromosome_path']:
            entity['disposition_reasons'].append('chromosome_extraction_absent')
        if not entity['whole_path']:
            entity['notes'].append('whole_assembly_file_absent')
        if entity['classifier_copy_path'] and not entity['classifier_copy_matches_chromosome']:
            entity['disposition_reasons'].append('classifier_reference_copy_differs_from_chromosome')
        if entity.get('legacy_exclude_marker'):
            entity['disposition_reasons'].append('legacy_exclude_marker_unresolved')
        if entity['crosswalk_review_flag'] == 'legacy_label_requires_review':
            entity['notes'].append('legacy_label_requires_review')
    return entities, conflicts


def public_duplicate_questions(entities, questions):
    """Identical chromosome sequence under different accessions: hold both, choose neither."""
    for slot, label in (('chromosome_sequence_sha256', 'chromosome'), ('whole_sequence_sha256', 'whole')):
        groups = defaultdict(list)
        for entity in entities.values():
            if entity[slot]:
                groups[entity[slot]].append(entity)
        for fingerprint, group in sorted(groups.items()):
            accessions = sorted({e['accession_version'] for e in group})
            if len(accessions) < 2:
                continue
            if len({e['assembly_numeric'] for e in group}) == 1:
                continue  # one assembly under two accession representations; see representation_questions
            question = questions.add(
                'public_discovery', f'identical_{label}_sequence_across_accessions',
                'blocking_study_cohort', accessions, len(group),
                f'identical sequence fingerprint {fingerprint[:16]}; '
                + '; '.join(f"{e['accession_version']} contigs={e[slot.replace('_sequence_sha256','_contigs')]} "
                            f"bp={e[slot.replace('_sequence_sha256','_total_bp')]}" for e in group),
                'results/handover-audit/2026-09-11/identical_sequence_groups.tsv',
                'NCBI assembly history: is one accession a replacement or re-release of the other, or are they '
                'independent submissions of the same material? Resolve in task 02 before the task 06 freeze. Do '
                'not prefer the larger or newer assembly without that history.')
            for entity in group:
                entity['provenance_question_ids'].append(question)
                entity['disposition_reasons'].append(f'duplicate_{label}_sequence_across_accessions')


def representation_questions(public, types, questions):
    """One assembly can appear as GCA and GCF, or as a study genome and a type anchor."""
    numerics = defaultdict(list)
    for entity in list(public.values()) + list(types.values()):
        if entity['assembly_numeric']:
            numerics[entity['assembly_numeric']].append(entity)
    for numeric, group in sorted(numerics.items()):
        accessions = sorted({e['accession_version'] for e in group})
        if len(accessions) > 1:
            question = questions.add(
                'public_discovery', 'multiple_accession_representations_of_one_assembly', 'informational',
                accessions, len(group),
                '; '.join(f"{e['accession_version']} [{e['cohort']}] {Path(e.get('path') or e.get('whole_path') or e.get('chromosome_path')).name}"
                          for e in group),
                'results/handover-audit/2026-09-11/assembly_manifest.tsv',
                'Carry one representation per assembly into the task 02 reference manifest and the task 06 '
                'freeze; record the GCA/GCF pairing and the version actually used rather than dropping one '
                'silently.')
            for entity in group:
                entity['provenance_question_ids'].append(question)
                entity['notes'].append('multiple_accession_representations_of_one_assembly')
        shared = sorted({e['entity_id'] for e in group if e['cohort'] == 'public_discovery'})
        if shared and any(e['cohort'] == 'type_reference_inventory' for e in group):
            for accession in shared:
                public[accession]['notes'].append('assembly_also_archived_as_a_type_reference')


def public_absence_questions(entities, questions):
    for reason, severity, resolution in (
            ('chromosome_extraction_absent', 'blocking_study_cohort',
             'Extract the chromosome from the archived whole assembly with the recorded method, or retrieve the '
             'assembly in task 02; do not substitute a different accession.'),
            ('no_local_sequence_for_roster_accession', 'blocking_study_cohort',
             'Retrieve the assembly in task 02 into a new directory, or record the roster row as withdrawn.'),
            ('legacy_exclude_marker_unresolved', 'blocking_naming',
             "Recover the reason behind the legacy 'Exclude?' marker from the analysis logbook before the entity "
             'is used or discarded.')):
        affected = sorted(e['entity_id'] for e in entities.values() if reason in e['disposition_reasons'])
        if affected:
            question = questions.add('public_discovery', reason, severity, affected, len(affected),
                                     f'{len(affected)} roster accessions', 'docs/audit/public-accession-crosswalk.tsv',
                                     resolution)
            for accession in affected:
                entities[accession]['provenance_question_ids'].append(question)
    missing_whole = sorted(e['entity_id'] for e in entities.values()
                           if 'whole_assembly_file_absent' in e['notes'])
    if missing_whole:
        question = questions.add(
            'public_discovery', 'whole_assembly_file_absent', 'informational', missing_whole,
            len(missing_whole), f'{len(missing_whole)} roster accessions have a chromosome extraction but no '
            'archived whole assembly', 'docs/audit/public-accession-crosswalk.tsv',
            'Retrieve the whole assemblies in task 02 if plasmid or whole-genome evidence is needed; the '
            'chromosome-based phylogeny is not blocked.')
        for accession in missing_whole:
            entities[accession]['provenance_question_ids'].append(question)
    review = sorted(e['entity_id'] for e in entities.values()
                    if 'legacy_label_requires_review' in e['notes'])
    if review:
        question = questions.add(
            'public_discovery', 'legacy_label_requires_review', 'blocking_naming', review, len(review),
            '; '.join(f"{a}={entities[a]['legacy_revised_species']}" for a in review),
            'docs/audit/public-accession-crosswalk.tsv',
            'Naming policy in task 12 after the type-anchored comparison; these entities stay in the study '
            'cohort but carry no accepted name.')
        for accession in review:
            entities[accession]['provenance_question_ids'].append(question)


def build_validation(audit_rows, questions, skani_evidence, accept_concordant):
    entities, unparsed = OrderedDict(), []
    grouped = defaultdict(lambda: defaultdict(list))
    for row in audit_rows:
        slot = VALIDATION_FILE_ROLES.get(row['role'])
        if not slot:
            continue
        parts = parse_validation_alias(row['filename'])
        if not parts:
            unparsed.append(dict(source=f'audit:{slot}', key=row['filename'], rows=1,
                                 conflict_kind='validation_filename_does_not_match_alias_grammar',
                                 disagreeing_fields='filename', values=row['filename'], detail=row['path']))
            continue
        grouped[parts['alias']][slot].append(dict(row, alias_parts=parts))

    for alias in sorted(grouped):
        slots = grouped[alias]
        candidates = sorted(slots.get('candidate', []), key=lambda row: row['path'])
        quarantine = sorted(slots.get('quarantine', []), key=lambda row: row['path'])
        misfiled = sorted(slots.get('misfiled_copy', []), key=lambda row: row['path'])
        parts = (candidates or quarantine or misfiled)[0]['alias_parts']
        concordant = len({row['sequence_multiset_sha256'] for row in candidates}) <= 1
        canonical = choose_canonical(candidates) if candidates else None
        entities[alias] = dict(
            entity_id=alias, cohort='clinical_validation', run_prefix=parts['run_prefix'],
            collection_date=parts['collection_date'], sample_number=parts['sample_number'],
            assembly_method=parts['assembly_method'], candidate_files=len(candidates),
            candidate_paths=join_field(candidates, 'path'), candidate_sha256=join_field(candidates, 'sha256'),
            candidate_sequence_sha256=join_field(candidates, 'sequence_multiset_sha256'),
            candidate_sequence_concordant=concordant if len(candidates) > 1 else '',
            canonical_path=canonical['path'] if canonical else '',
            canonical_sha256=canonical['sha256'] if canonical else '',
            canonical_sequence_sha256=canonical['sequence_multiset_sha256'] if canonical else '',
            canonical_contigs=canonical['contigs'] if canonical else '',
            canonical_total_bp=canonical['total_bp'] if canonical else '',
            canonical_qc_flags=canonical['qc_flags'] if canonical else '',
            quarantine_files=len(quarantine), quarantine_paths=join_field(quarantine, 'path'),
            quarantine_sequence_matches_candidate=all(
                row['sequence_multiset_sha256'] in {c['sequence_multiset_sha256'] for c in candidates}
                for row in quarantine) if quarantine and candidates else '',
            misfiled_copy_paths=join_field(misfiled, 'path'),
            files=len(candidates) + len(quarantine) + len(misfiled),
            entity_role='validation_candidate', usage_restrictions=[],
            selected_inference_entity_id='', alternate_entity_ids=[], adjudication_ids=[],
            disposition='', disposition_reasons=[], notes=[], provenance_question_ids=[])

    for alias, entity in entities.items():
        slots = grouped[alias]
        candidates = sorted(slots.get('candidate', []), key=lambda row: row['path'])
        quarantine = sorted(slots.get('quarantine', []), key=lambda row: row['path'])
        misfiled = sorted(slots.get('misfiled_copy', []), key=lambda row: row['path'])
        if not candidates:
            entity['disposition_reasons'].append('alias_present_only_outside_the_validation_set')
        if len(candidates) > 1:
            concordant = entity['candidate_sequence_concordant'] is True
            rule = ('alias_collision_sequence_concordant' if concordant
                    else 'alias_collision_sequence_discordant')
            evidence = ('; '.join(f"{row['filename']} contigs={row['contigs']} bp={row['total_bp']} "
                                  f"seq={row['sequence_multiset_sha256'][:12]}" for row in candidates))
            reuse = pair_evidence(skani_evidence, candidates)
            question = questions.add(
                'clinical_validation', rule,
                'blocking_validation_denominator', [alias], len(candidates),
                evidence + ('; ' + reuse if reuse else ''),
                'results/handover-audit/2026-09-11/validation_collisions.tsv'
                + ('; reused skani/minimap2 comparison' if reuse else ''),
                "Confirm against the clinical sample sheet which filename is the isolate's record."
                + ('' if concordant else
                   ' The two files are not the same assembly: decide which assembly run belongs to this isolate '
                   'before any accuracy claim uses it.'))
            entity['provenance_question_ids'].append(question)
            if concordant and accept_concordant:
                entity['notes'].append('alias_collision_sequence_concordant_canonical_file_selected')
            else:
                entity['disposition_reasons'].append(rule)
        if quarantine and entity['quarantine_sequence_matches_candidate'] is False:
            differing = [row for row in quarantine
                         if row['sequence_multiset_sha256'] not in
                         {c['sequence_multiset_sha256'] for c in candidates}]
            question = questions.add(
                'clinical_validation', 'quarantined_copy_sequence_differs_from_retained_candidate',
                'blocking_validation_denominator', [alias], len(differing),
                '; '.join(f"{row['filename']} contigs={row['contigs']} bp={row['total_bp']} "
                          f"seq={row['sequence_multiset_sha256'][:12]}" for row in differing)
                + f"; retained candidate seq={entity['canonical_sequence_sha256'][:12]}",
                'results/handover-audit/2026-09-11/assembly_manifest.tsv',
                'The archived quarantine was treated as a leading-zero duplicate, but this file is a different '
                'assembly. Establish which assembly the clinical record refers to before the alias is used.')
            entity['provenance_question_ids'].append(question)
            entity['disposition_reasons'].append('quarantined_copy_sequence_differs_from_retained_candidate')
        elif quarantine:
            entity['notes'].append('quarantined_leading_zero_copy_sequence_matches_candidate')
        if misfiled:
            question = questions.add(
                'clinical_validation', 'clinical_query_stored_in_public_chromosome_directory',
                'informational', [alias], len(misfiled),
                '; '.join(row['path'] for row in misfiled),
                'results/handover-audit/2026-09-11/assembly_manifest.tsv',
                'Leave the archived file in place and address it with explicit manifests: a wildcard over the '
                'public chromosome directory would ingest a clinical query as a reference genome.')
            entity['provenance_question_ids'].append(question)
            entity['notes'].append('clinical_query_copy_present_in_public_chromosome_directory')

    cross = defaultdict(list)
    for entity in entities.values():
        for fingerprint in {h for h in entity['candidate_sequence_sha256'].split(';') if h}:
            cross[fingerprint].append(entity)
    for fingerprint, group in sorted(cross.items()):
        aliases = sorted({e['entity_id'] for e in group})
        if len(aliases) < 2:
            continue
        byte_identical = len({e['canonical_sha256'] for e in group}) == 1
        question = questions.add(
            'clinical_validation',
            'identical_sequence_under_different_sample_aliases', 'blocking_validation_denominator',
            aliases, sum(e['candidate_files'] for e in group),
            f'identical sequence fingerprint {fingerprint[:16]}'
            + ('; the files are also byte-identical' if byte_identical else ''),
            'results/handover-audit/2026-09-11/identical_sequence_groups.tsv',
            'Different sample numbers carry the same assembly. Confirm from the sample sheet whether one number '
            'is a transcription error or a file was copied over another; do not treat the aliases as one isolate '
            'and do not assume both isolates exist.')
        for entity in group:
            entity['provenance_question_ids'].append(question)
            entity['disposition_reasons'].append('identical_sequence_under_different_sample_aliases')
    return entities, unparsed


def build_type_inventory(audit_rows, questions):
    entities = OrderedDict()
    for row in sorted((r for r in audit_rows if r['role'] == 'type_reference_candidate'),
                      key=lambda r: r['path']):
        meta = parse_accession(row['accession_version']) or dict(
            accession_db='', assembly_numeric='', accession_version='')
        entities[row['path']] = dict(
            entity_id=row['path'], cohort='type_reference_inventory',
            accession_version=row['accession_version'], accession_db=meta['accession_db'],
            assembly_numeric=meta['assembly_numeric'], filename=row['filename'], path=row['path'],
            sha256=row['sha256'], sequence_sha256=row['sequence_multiset_sha256'],
            contigs=row['contigs'], total_bp=row['total_bp'], qc_flags=row['qc_flags'],
            first_header=row['first_header'], files=1,
            entity_role='type_reference_inventory_record', usage_restrictions=[],
            selected_inference_entity_id='', alternate_entity_ids=[], adjudication_ids=[],
            disposition='defer', disposition_reasons=['type_panel_verified_in_task_02'],
            notes=['species_named_file_without_accession'] if not row['accession_version'] else [],
            provenance_question_ids=[])
    return entities


def load_adjudications(paths):
    """Read and validate accepted adjudication tables; never guess at an unreadable decision."""
    decisions = []
    for path in paths:
        rows, fields = read_tsv_with_fields(path)
        require_fields(fields, ADJUDICATION_FIELDS, 'adjudication table', path)
        require_unique(rows, 'decision_id', 'adjudication table', path)
        decisions.extend(dict(row, decision_source=str(path)) for row in rows)
    require_unique(decisions, 'decision_id', 'adjudication tables',
                   ' + '.join(str(path) for path in paths))
    for decision in decisions:
        actions = [action.strip() for action in decision['proposed_cohort_action'].split(';')
                   if action.strip()]
        unknown = sorted(set(actions) - KNOWN_ACTIONS)
        if unknown:
            raise ReconcileError(
                f"{decision['decision_id']}: unknown cohort action {', '.join(unknown)}; refusing "
                'to apply a decision this script does not implement')
        if not actions:
            raise ReconcileError(f"{decision['decision_id']}: no cohort action recorded")
        review = decision['coordinator_review'].strip()
        if review.startswith(REVIEW_ACCEPTED):
            state = 'accepted'
        elif review.startswith(REVIEW_NO_ADJUDICATION):
            state = 'no_adjudication_required'
        elif review.startswith(REVIEW_WITHHELD):
            state = 'withheld'
        else:
            raise ReconcileError(f"{decision['decision_id']}: unrecognised coordinator review "
                                 f"state '{review}'")
        decision.update(actions=actions, review_state=state,
                        applied=state in {'accepted', 'no_adjudication_required'},
                        questions=expand_question_ids(decision['question_ids']))
    return sorted(decisions, key=lambda decision: decision['decision_id'])


def decision_entities(decision, public):
    """Exact accession.version entities that must already exist in the public cohort."""
    listed = [entity.strip() for entity in decision['entities'].split(';') if entity.strip()]
    accessions = [entity for entity in listed if ACCESSION.match(entity)]
    if len(accessions) != len(listed):
        if set(decision['actions']) <= REPRESENTATION_ACTIONS:
            return []  # a cohort-wide representation policy, not an entity list
        raise ReconcileError(f"{decision['decision_id']}: entity list is not accession.version: "
                             f"{decision['entities']}")
    unknown = [entity for entity in accessions if entity not in public]
    if unknown:
        raise ReconcileError(f"{decision['decision_id']}: names entities that are not public "
                             f"discovery accessions: {', '.join(unknown)}")
    return accessions


def resolve_canonical(decision, entity_ids, public):
    """Map the accepted canonical database record onto exactly one cohort entity."""
    canonical = decision['canonical_record_if_accepted'].strip()
    record = parse_accession(canonical)
    if not record:
        raise ReconcileError(f"{decision['decision_id']}: canonical record '{canonical}' is not an "
                             'accession.version')
    matches = [entity for entity in entity_ids
               for other in [parse_accession(entity)]
               if other and other['assembly_numeric'] == record['assembly_numeric']
               and other['accession_version_number'] == record['accession_version_number']]
    if len(matches) != 1:
        raise ReconcileError(
            f"{decision['decision_id']}: canonical record {canonical} matches {len(matches)} of "
            f"the listed entities ({', '.join(entity_ids)}); an ambiguous canonical record is not "
            'an adjudication')
    return matches[0], canonical


def apply_adjudications(decisions, public, questions):
    """Apply validated decisions to public dispositions and annotate the question ledger."""
    applied = []
    for decision in decisions:
        entity_ids = decision_entities(decision, public)
        kind, detail, scope = 'recorded_only', '', set(entity_ids) or None
        if not decision['applied']:
            questions.annotate(decision['questions'], 'deferred', decision['decision_id'],
                               f"coordinator review {decision['coordinator_review']}", scope)
        elif set(decision['actions']) <= PAIR_ACTIONS:
            canonical_entity, canonical_record = resolve_canonical(decision, entity_ids, public)
            alternates = [entity for entity in entity_ids if entity != canonical_entity]
            if not alternates:
                raise ReconcileError(f"{decision['decision_id']}: a duplicate-submission decision "
                                     'must name at least one alternate submission')
            selected = public[canonical_entity]
            selected['disposition'] = 'keep'
            selected['disposition_reasons'] = [
                reason for reason in selected['disposition_reasons']
                if not reason.startswith('duplicate_')]
            selected['notes'].append('selected_inference_entity_for_duplicate_submission_pair')
            selected['canonical_public_record'] = canonical_record
            selected['alternate_entity_ids'].extend(alternates)
            selected['adjudication_ids'].append(decision['decision_id'])
            for accession in alternates:
                alternate = public[accession]
                alternate['disposition'] = 'exclude'
                alternate['disposition_reasons'].append(
                    'duplicate_submission_of_selected_inference_entity')
                alternate['selected_inference_entity_id'] = canonical_entity
                alternate['canonical_public_record'] = canonical_record
                alternate['entity_role'] = 'alternate_submission_record'
                alternate['usage_restrictions'].append('not_an_independent_inference_unit')
                alternate['notes'].append('retained_in_provenance_and_exclusions')
                alternate['adjudication_ids'].append(decision['decision_id'])
            kind = 'one_inference_unit_per_duplicate_submission_pair'
            detail = (f'{canonical_entity} selected as the inference entity for '
                      f'{canonical_record}; alternates ' + ';'.join(alternates))
            questions.annotate(decision['questions'], 'accepted', decision['decision_id'],
                               detail, scope)
        elif set(decision['actions']) <= EXCLUDE_ACTIONS:
            for accession in entity_ids:
                entity = public[accession]
                entity['disposition'] = 'exclude'
                entity['disposition_reasons'].append(
                    f"adjudicated_exclusion:{decision['proposed_relationship']}")
                entity['entity_role'] = 'excluded_record'
                entity['usage_restrictions'].extend(
                    ['excluded_from_phylogeny', 'excluded_from_classifier_reference'])
                entity['notes'].append('audit_and_provenance_record_retained')
                entity['adjudication_ids'].append(decision['decision_id'])
            kind = 'excluded_from_phylogeny_and_classifier_reference'
            detail = ', '.join(entity_ids)
            questions.annotate(decision['questions'], 'accepted', decision['decision_id'],
                               detail, scope)
        else:
            crosswalked = sorted(
                entity['entity_id'] for entity in public.values()
                if 'multiple_accession_representations_of_one_assembly' in entity['notes'])
            for accession in crosswalked:
                entity = public[accession]
                entity['notes'].append('one_assembly_inference_unit_accession_crosswalk_retained')
                entity['adjudication_ids'].append(decision['decision_id'])
            kind = 'one_assembly_inference_unit_with_retained_crosswalk'
            detail = (f'{len(crosswalked)} public entities carry a retained accession crosswalk; '
                      'the paired representation adds no inference unit')
            questions.annotate(decision['questions'], 'resolved', decision['decision_id'], detail)
        applied.append(dict(decision_id=decision['decision_id'],
                            review_state=decision['review_state'],
                            applied=decision['applied'], action_kind=kind,
                            entities=';'.join(entity_ids), questions=';'.join(decision['questions']),
                            canonical_record=decision['canonical_record_if_accepted'],
                            detail=detail, decision_source=decision['decision_source']))
    return applied


def load_resolved_files(paths, public):
    """Read and validate resolved file-gap manifests against the public cohort."""
    rows = []
    for path in paths:
        loaded, fields = read_tsv_with_fields(path)
        require_fields(fields, RESOLVED_FILE_FIELDS, 'resolved file manifest', path)
        require_unique(loaded, 'accession', 'resolved file manifest', path)
        rows.extend(dict(row, resolved_source=str(path)) for row in loaded)
    require_unique(rows, 'accession', 'resolved file manifests',
                   ' + '.join(str(path) for path in paths))
    for row in rows:
        if not ACCESSION.match(row['accession'].strip()):
            raise ReconcileError(f"resolved file manifest key '{row['accession']}' is not an "
                                 'accession.version')
        if row['accession'] not in public:
            raise ReconcileError(f"resolved file manifest names {row['accession']}, which is not a "
                                 'public discovery entity')
        if row['validation'].strip() != 'pass':
            raise ReconcileError(f"{row['accession']}: resolved file record is marked "
                                 f"'{row['validation']}', not 'pass'")
    return sorted(rows, key=lambda row: row['accession'])


def verify_resolved_file_row(row):
    """Re-verify paths, file checksums, sequence hash, length and named chromosome accession."""
    checks = []
    for label, path_field, hash_field, origin_field in (
            ('whole', 'whole_path', 'whole_sha256', 'whole_origin'),
            ('chromosome', 'chromosome_path', 'chromosome_file_sha256', 'chromosome_origin')):
        path = row[path_field].strip()
        origin = row[origin_field].strip()
        if not path:
            if origin:
                raise ReconcileError(f"{row['accession']}: {label} origin '{origin}' without a path")
            continue
        if origin not in ARCHIVE_ORIGINS and not origin.startswith(RETRIEVED_ORIGIN_PREFIXES):
            raise ReconcileError(f"{row['accession']}: unknown {label} origin '{origin}'")
        target = Path(path)
        if not target.is_file():
            raise ReconcileError(f"{row['accession']}: {label} file is missing: {path}")
        observed = digest_file(target)
        if observed != row[hash_field].strip():
            raise ReconcileError(f"{row['accession']}: {label} file checksum {observed} does not "
                                 f"match the recorded {row[hash_field].strip()}")
        checks.append(f'{label}:{origin}:file_sha256_verified')
    chromosome = row['chromosome_path'].strip()
    if chromosome:
        records = read_fasta_records(chromosome)
        if len(records) != 1:
            raise ReconcileError(f"{row['accession']}: a named chromosome file must hold exactly "
                                 f'one record, found {len(records)}')
        header, sequence = records[0]
        observed = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        if observed != row['chromosome_sequence_sha256'].strip():
            raise ReconcileError(f"{row['accession']}: chromosome sequence checksum {observed} does "
                                 f"not match the recorded {row['chromosome_sequence_sha256'].strip()}")
        if str(len(sequence)) != row['chromosome_length'].strip():
            raise ReconcileError(f"{row['accession']}: chromosome length {len(sequence)} does not "
                                 f"match the recorded {row['chromosome_length'].strip()}")
        if header.split()[0] != row['chromosome_accession'].strip():
            raise ReconcileError(f"{row['accession']}: chromosome file names "
                                 f"{header.split()[0]}, not {row['chromosome_accession'].strip()}")
        checks.append('chromosome:sequence_sha256_length_and_named_accession_verified')
        row['chromosome_multiset_sha256'] = sequence_multiset_sha256(records)
    return checks


def resolved_file_conflicts(row, entity):
    """A resolved record may complete the audit; it may never quietly contradict or replace it."""
    conflicts = []
    for label, audit_path, audit_hash, path_field, hash_field, origin_field in (
            ('chromosome', entity['chromosome_path'], entity['chromosome_sha256'],
             'chromosome_path', 'chromosome_file_sha256', 'chromosome_origin'),
            ('whole', entity['whole_path'], entity['whole_sha256'],
             'whole_path', 'whole_sha256', 'whole_origin')):
        path, origin = row[path_field].strip(), row[origin_field].strip()
        if not path:
            continue
        if origin in ARCHIVE_ORIGINS:
            if not audit_path:
                conflicts.append(f'{label}: an archived file the FASTA audit did not record')
            elif audit_path != path:
                conflicts.append(f'{label}: archived path disagrees with the audit '
                                 f'({audit_path} vs {path})')
            elif audit_hash != row[hash_field].strip():
                conflicts.append(f'{label}: archived checksum disagrees with the audit')
        elif audit_path:
            conflicts.append(f'{label}: a retrieved or derived file would compete with the '
                             f'archived file the audit already recorded ({audit_path})')
    return conflicts


def apply_resolved_files(rows, public, questions, source_label):
    """Attach verified file-gap representations without folding them into the audited roles."""
    applied = []
    for row in rows:
        entity = public[row['accession']]
        checks = verify_resolved_file_row(row)
        conflicts = resolved_file_conflicts(row, entity)
        if conflicts:
            raise ReconcileError(f"{row['accession']}: resolved file record conflicts with the "
                                 'FASTA audit: ' + '; '.join(conflicts))
        archived_before = bool(entity['chromosome_path'] or entity['whole_path'])
        entity.update(
            resolved_whole_path=row['whole_path'].strip(),
            resolved_whole_origin=row['whole_origin'].strip(),
            resolved_whole_sha256=row['whole_sha256'].strip(),
            resolved_chromosome_path=row['chromosome_path'].strip(),
            resolved_chromosome_origin=row['chromosome_origin'].strip(),
            resolved_chromosome_file_sha256=row['chromosome_file_sha256'].strip(),
            resolved_chromosome_sequence_sha256=row['chromosome_sequence_sha256'].strip(),
            resolved_chromosome_multiset_sha256=row.get('chromosome_multiset_sha256', ''),
            resolved_chromosome_accession=row['chromosome_accession'].strip(),
            resolved_chromosome_length=row['chromosome_length'].strip(),
            resolved_chromosome_header=row['chromosome_header'].strip(),
            resolved_ncbi_genome_notes=row['ncbi_genome_notes'].strip(),
            resolved_file_checks=';'.join(checks),
            resolved_file_source=row['resolved_source'])
        for reason in ('chromosome_extraction_absent', 'no_local_sequence_for_roster_accession'):
            if reason in entity['disposition_reasons']:
                entity['disposition_reasons'].remove(reason)
                entity['notes'].append(f'{reason}_closed_by_verified_file_record')
        if 'whole_assembly_file_absent' in entity['notes'] and entity['resolved_whole_path']:
            entity['notes'].append('whole_assembly_representation_verified_from_resolved_record')
        if not archived_before and UNRESOLVED_TAXONOMY_NOTE in entity['resolved_ncbi_genome_notes']:
            # No archived sequence of its own and an unresolved source organism: keep the sequence
            # for placement testing, but never as a named reference before tasks 04/05/07-10.
            entity['entity_role'] = 'context_candidate'
            entity['usage_restrictions'].extend(['no_named_classifier_reference',
                                                 'no_named_reference_accuracy_denominator',
                                                 'placement_pending_tasks_04_05_07_10'])
            entity['notes'].append('unresolved_taxonomy_context_candidate')
        applied.append(dict(accession=row['accession'], entity_role=entity['entity_role'],
                            whole_origin=entity['resolved_whole_origin'],
                            chromosome_origin=entity['resolved_chromosome_origin'],
                            chromosome_accession=entity['resolved_chromosome_accession'],
                            checks=entity['resolved_file_checks'],
                            ncbi_genome_notes=entity['resolved_ncbi_genome_notes'],
                            questions=';'.join(expand_question_ids(row['question_ids'])),
                            resolved_source=row['resolved_source']))
    covered = {row['accession'] for row in rows}
    if covered:
        questions.annotate(sorted({question for row in rows
                                   for question in expand_question_ids(row['question_ids'])}),
                           'resolved', source_label,
                           'verified local file representation with recorded origin and checksums',
                           covered)
    return applied


def manifest_changes(before_rows, after_rows, manifest, key='entity_id',
                     compare=('disposition', 'disposition_reasons', 'entity_role')):
    """Row-level before/after diff for one shared manifest."""
    before = {row[key]: row for row in before_rows}
    after = {row[key]: row for row in after_rows}
    changes, retained = [], 0
    # A column the previous manifest did not have is a schema change, not 1,059 row changes.
    if before_rows and after_rows:
        changes.extend(dict(manifest=manifest, entity_id='', change='column_added', field=field,
                            before='', after='')
                       for field in sorted(set(after_rows[0]) - set(before_rows[0])))
        changes.extend(dict(manifest=manifest, entity_id='', change='column_removed', field=field,
                            before='', after='')
                       for field in sorted(set(before_rows[0]) - set(after_rows[0])))
        compare = tuple(field for field in compare
                        if field in set(before_rows[0]) & set(after_rows[0]))
    for entity_id in sorted(set(before) | set(after)):
        if entity_id not in before:
            changes.append(dict(manifest=manifest, entity_id=entity_id, change='added', field='',
                                before='', after=after[entity_id].get('disposition', '')))
        elif entity_id not in after:
            changes.append(dict(manifest=manifest, entity_id=entity_id, change='removed', field='',
                                before=before[entity_id].get('disposition', ''), after=''))
        else:
            differing = [field for field in compare
                         if str(before[entity_id].get(field, '')) != str(after[entity_id].get(field, ''))]
            if differing:
                changes.extend(dict(manifest=manifest, entity_id=entity_id, change='changed',
                                    field=field, before=before[entity_id].get(field, ''),
                                    after=after[entity_id].get(field, '')) for field in differing)
            else:
                retained += 1
    return changes, retained


def finalize(entities, keep_label='keep'):
    for entity in entities.values():
        if not entity['disposition']:
            entity['disposition'] = 'hold' if entity['disposition_reasons'] else keep_label
        for field in ('disposition_reasons', 'notes', 'provenance_question_ids',
                      'usage_restrictions', 'alternate_entity_ids', 'adjudication_ids'):
            entity[field] = ';'.join(dict.fromkeys(entity[field])) if isinstance(entity[field], list) \
                else entity[field]
    return entities


def transitions(stage, rule, unit, count, note):
    return dict(stage=stage, rule=rule, unit=unit, count=count, note=note)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--audit', type=Path, required=True, help='assembly_manifest.tsv from the FASTA audit')
    parser.add_argument('--crosswalk', type=Path, required=True, help='public accession crosswalk (roster)')
    parser.add_argument('--legacy-metadata', type=Path, required=True, help='legacy reference_metadata.tsv')
    parser.add_argument('--output', type=Path, required=True, help='fresh results/cohort/RUN_ID directory')
    parser.add_argument('--manifest-dir', type=Path, default=Path('manifests'))
    parser.add_argument('--skani-evidence', type=Path, action='append', default=[],
                        help='skani dist table reused as pair evidence; may be repeated')
    parser.add_argument('--accept-concordant-alias-collisions', action='store_true',
                        help='keep aliases whose colliding files carry the same sequence, recording the '
                             'canonical filename choice; off by default so the denominator stays conservative')
    parser.add_argument('--adjudications', type=Path, action='append', default=[],
                        help='accepted provenance adjudication table; may be repeated')
    parser.add_argument('--resolved-files', type=Path, action='append', default=[],
                        help='resolved file-gap manifest of verified local representations; '
                             'may be repeated')
    parser.add_argument('--defer-clinical-questions', metavar='SOURCE',
                        help='annotate every still-open clinical question as deferred with this '
                             'decision source; dispositions, counts and manifests are unchanged')
    parser.add_argument('--stage-only', action='store_true',
                        help='write the proposed manifests into the run directory and leave the '
                             'shared manifest directory untouched')
    parser.add_argument('--force', action='store_true', help='allow a non-empty output directory')
    args = parser.parse_args(argv)

    output = args.output.resolve()
    if output.exists() and any(output.iterdir()) and not args.force:
        parser.error('Use a fresh output directory; pass --force only to rewrite this run')
    output.mkdir(parents=True, exist_ok=True)

    audit_rows = read_tsv(args.audit)
    crosswalk = read_tsv(args.crosswalk)
    legacy = read_tsv(args.legacy_metadata)
    skani_evidence = load_skani_evidence(args.skani_evidence)

    questions = Questions()
    public, conflicts = build_public(audit_rows, crosswalk, legacy, questions)
    types = build_type_inventory(audit_rows, questions)
    public_duplicate_questions(public, questions)
    representation_questions(public, types, questions)
    public_absence_questions(public, questions)
    validation, unparsed = build_validation(audit_rows, questions, skani_evidence,
                                            args.accept_concordant_alias_collisions)
    conflicts.extend(unparsed)

    # Validate every optional decision and file record before any of them is applied.
    decisions = load_adjudications(args.adjudications)
    resolved_rows = load_resolved_files(args.resolved_files, public)
    applied_files = apply_resolved_files(
        resolved_rows, public, questions,
        ' + '.join(str(path) for path in args.resolved_files) or 'resolved-files')
    applied_decisions = apply_adjudications(decisions, public, questions)
    if args.defer_clinical_questions:
        questions.annotate([row['question_id'] for row in questions.rows
                            if row['cohort'] == 'clinical_validation' and row['status'] == 'open'],
                           'deferred', args.defer_clinical_questions,
                           'clinical validation cohort parked; dispositions and counts unchanged')
    finalize(public)
    finalize(validation)
    finalize(types)

    role_counts = Counter(row['role'] for row in audit_rows)
    assigned_files = (sum(e['files'] for e in public.values())
                      + sum(e['files'] for e in validation.values())
                      + sum(e['files'] for e in types.values()))
    if assigned_files != len(audit_rows):
        raise ValueError(f'{len(audit_rows) - assigned_files} audited files were not assigned to an entity')

    rows = [transitions('input', 'audited_assembly_files', 'files', len(audit_rows), str(args.audit)),
            transitions('input', 'discovery_roster_rows', 'accessions', len(crosswalk), str(args.crosswalk)),
            transitions('input', 'legacy_metadata_rows', 'genomes', len(legacy), str(args.legacy_metadata))]
    for role, count in sorted(role_counts.items()):
        rows.append(transitions('input', f'audit_role:{role}', 'files', count, 'role assigned by the FASTA audit'))
    for cohort, entities, unit in (('public_discovery', public, 'accessions'),
                                   ('clinical_validation', validation, 'sample_aliases'),
                                   ('type_reference_inventory', types, 'files')):
        rows.append(transitions('entities', f'{cohort}:entities', unit, len(entities),
                                'files collapsed onto stable identities'))
        rows.append(transitions('entities', f'{cohort}:files', 'files',
                                sum(e['files'] for e in entities.values()), 'every audited file is assigned once'))
        reasons = Counter(reason for e in entities.values()
                          for reason in e['disposition_reasons'].split(';') if reason)
        for reason, count in sorted(reasons.items()):
            rows.append(transitions('rules', f'{cohort}:reason:{reason}', unit, count,
                                    'counted per entity; an entity may carry several reasons'))
        notes = Counter(note for e in entities.values() for note in e['notes'].split(';') if note)
        for note, count in sorted(notes.items()):
            rows.append(transitions('rules', f'{cohort}:note:{note}', unit, count, 'recorded, not disqualifying'))
        for disposition, count in sorted(Counter(e['disposition'] for e in entities.values()).items()):
            rows.append(transitions('final', f'{cohort}:{disposition}', unit, count, 'disposition after all rules'))
    rows.append(transitions('adjudications', 'decisions_read', 'decisions', len(decisions),
                            ' + '.join(str(path) for path in args.adjudications) or 'none supplied'))
    rows.append(transitions('adjudications', 'decisions_applied', 'decisions',
                            sum(1 for decision in applied_decisions if decision['applied']),
                            'accepted or verified-bookkeeping decisions actually applied'))
    for decision in applied_decisions:
        rows.append(transitions('adjudications', f"decision:{decision['decision_id']}", 'entities',
                                len([e for e in decision['entities'].split(';') if e]),
                                f"{decision['review_state']}; {decision['action_kind']}"))
    rows.append(transitions('resolved_files', 'records_verified', 'accessions', len(applied_files),
                            'paths, file checksums, sequence hash, length and named accession '
                            're-checked on disk'))
    for field in ('whole_origin', 'chromosome_origin'):
        for origin, count in sorted(Counter(record[field] for record in applied_files
                                            if record[field]).items()):
            rows.append(transitions('resolved_files', f'{field}:{origin}', 'accessions', count,
                                    'additional sequence representation, not an audited file'))
    for role, count in sorted(Counter(record['entity_role'] for record in applied_files).items()):
        rows.append(transitions('resolved_files', f'entity_role:{role}', 'accessions', count,
                                'role assigned from the verified file record'))
    for status, count in sorted(Counter(row['status'] for row in questions.rows).items()):
        rows.append(transitions('questions', f'status:{status}', 'questions', count,
                                'every question is retained; status is annotated, never erased'))
    rows.append(transitions('final', 'key_conflicts', 'conflicts', len(conflicts),
                            'repeated keys reported rather than silently collapsed'))
    rows.append(transitions('final', 'provenance_questions', 'questions', len(questions.rows),
                            'open identity questions carried to the coordinator'))
    rows.append(transitions('final', 'study_candidates', 'accessions',
                            sum(e['disposition'] == 'keep' for e in public.values()), 'manifests/study-candidates.tsv'))
    rows.append(transitions('final', 'validation_candidates', 'sample_aliases',
                            sum(e['disposition'] == 'keep' for e in validation.values()),
                            'manifests/validation-candidates.tsv'))

    public_fields = [f for f in next(iter(public.values()))]
    validation_fields = [f for f in next(iter(validation.values()))]
    type_fields = [f for f in next(iter(types.values()))]
    write_tsv(output / 'public-entities.tsv', public.values(), public_fields)
    write_tsv(output / 'validation-entities.tsv', validation.values(), validation_fields)
    write_tsv(output / 'type-reference-entities.tsv', types.values(), type_fields)
    write_tsv(output / 'provenance-questions.tsv', questions.rows,
              ['question_id', 'cohort', 'rule', 'severity', 'entities', 'entity_ids', 'file_count',
               'evidence', 'evidence_source', 'resolution_required', 'status', 'status_source',
               'status_detail'])
    write_tsv(output / 'conflicts.tsv', conflicts,
              ['source', 'key', 'rows', 'conflict_kind', 'disagreeing_fields', 'values', 'detail'])
    write_tsv(output / 'applied-adjudications.tsv', applied_decisions,
              ['decision_id', 'review_state', 'applied', 'action_kind', 'entities', 'questions',
               'canonical_record', 'detail', 'decision_source'])
    write_tsv(output / 'applied-resolved-files.tsv', applied_files,
              ['accession', 'entity_role', 'whole_origin', 'chromosome_origin',
               'chromosome_accession', 'checks', 'ncbi_genome_notes', 'questions',
               'resolved_source'])

    exclusion_fields = ['entity_id', 'cohort', 'disposition', 'disposition_reasons', 'entity_role',
                        'usage_restrictions', 'selected_inference_entity_id',
                        'canonical_public_record', 'adjudication_ids', 'notes',
                        'provenance_question_ids', 'files', 'representative_path', 'sha256']
    exclusions = []
    for entities, path_field, hash_field in ((public, 'chromosome_path', 'chromosome_sha256'),
                                             (validation, 'canonical_path', 'canonical_sha256'),
                                             (types, 'path', 'sha256')):
        exclusions.extend(dict(e, representative_path=e.get(path_field, '') or e.get('whole_path', '')
                               or e.get('resolved_chromosome_path', ''),
                               sha256=e.get(hash_field, '') or e.get('whole_sha256', '')
                               or e.get('resolved_chromosome_file_sha256', ''))
                          for e in entities.values() if e['disposition'] != 'keep')
    proposed = OrderedDict((
        ('study-candidates.tsv',
         ([e for e in public.values() if e['disposition'] == 'keep'], public_fields)),
        ('validation-candidates.tsv',
         ([e for e in validation.values() if e['disposition'] == 'keep'], validation_fields)),
        ('exclusions.tsv', (exclusions, exclusion_fields))))

    # Stage and diff first; the shared manifests are replaced only after the proposed rows exist.
    manifests = args.manifest_dir.resolve()
    staged = output / 'proposed-manifests'
    diff_rows, diff_counts = [], {}
    for name, (entries, fields) in proposed.items():
        before = read_tsv(manifests / name) if (manifests / name).exists() else []
        write_tsv(staged / name, entries, fields)
        after = read_tsv(staged / name)
        changes, retained = manifest_changes(before, after, name)
        diff_rows.extend(changes)
        counted = dict(Counter(change['change'] for change in changes), retained=retained,
                       before=len(before), after=len(after))
        diff_counts[name] = counted
        for change, count in sorted(counted.items()):
            rows.append(transitions('manifest_diff', f'{name}:{change}', 'rows', count,
                                    'shared manifest before this run versus the proposed rows'))
    write_tsv(output / 'manifest-diff.tsv', diff_rows,
              ['manifest', 'entity_id', 'change', 'field', 'before', 'after'])
    if not args.stage_only:
        for name, (entries, fields) in proposed.items():
            write_tsv(manifests / name, entries, fields)

    write_tsv(output / 'count-transitions.tsv', rows, ['stage', 'rule', 'unit', 'count', 'note'])

    produced = ['public-entities.tsv', 'validation-entities.tsv', 'type-reference-entities.tsv',
                'count-transitions.tsv', 'provenance-questions.tsv', 'conflicts.tsv',
                'applied-adjudications.tsv', 'applied-resolved-files.tsv', 'manifest-diff.tsv']
    produced.extend(f'proposed-manifests/{name}' for name in MANIFEST_NAMES)
    summary = dict(
        created_utc=datetime.now(timezone.utc).isoformat(), command=sys.argv,
        script_sha256=digest_file(__file__),
        inputs={str(p): digest_file(p) for p in [args.audit, args.crosswalk, args.legacy_metadata]
                + list(args.skani_evidence) + list(args.adjudications) + list(args.resolved_files)},
        accept_concordant_alias_collisions=args.accept_concordant_alias_collisions,
        defer_clinical_questions=args.defer_clinical_questions or '',
        stage_only=args.stage_only, shared_manifests_written=not args.stage_only,
        adjudications={'read': len(decisions),
                       'applied': sum(1 for d in applied_decisions if d['applied']),
                       'by_review_state': dict(Counter(d['review_state'] for d in applied_decisions)),
                       'by_action_kind': dict(Counter(d['action_kind'] for d in applied_decisions))},
        resolved_file_records={'verified': len(applied_files),
                               'by_chromosome_origin': dict(Counter(
                                   r['chromosome_origin'] for r in applied_files if r['chromosome_origin'])),
                               'by_whole_origin': dict(Counter(
                                   r['whole_origin'] for r in applied_files if r['whole_origin']))},
        manifest_diff=diff_counts,
        audited_files=len(audit_rows), assigned_files=assigned_files,
        entities={'public_discovery': len(public), 'clinical_validation': len(validation),
                  'type_reference_inventory': len(types)},
        dispositions={cohort: dict(Counter(e['disposition'] for e in entities.values()))
                      for cohort, entities in (('public_discovery', public),
                                               ('clinical_validation', validation),
                                               ('type_reference_inventory', types))},
        provenance_questions=len(questions.rows),
        provenance_questions_by_severity=dict(Counter(q['severity'] for q in questions.rows)),
        provenance_questions_by_status=dict(Counter(q['status'] for q in questions.rows)),
        key_conflicts=len(conflicts), caveats=CAVEATS)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (output / 'outputs.sha256').write_text(''.join(
        f'{digest_file(output / name)}  {name}\n' for name in produced)
        + ('' if args.stage_only else
           ''.join(f'{digest_file(manifests / name)}  {manifests / name}\n'
                   for name in MANIFEST_NAMES)))
    return summary


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
