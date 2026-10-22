#!/usr/bin/env python3
"""Bounded descriptive Standard Kraken comparison; never infers genome hard calls."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import time
from pathlib import Path

from apply_reporting_categories import make_crosswalk, parse_label

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'CONTEXT.md').is_file() and (p / 'scripts').is_dir())
MANIFEST = ROOT / 'results/benchmark-design/2026-10-03/task14-astra-whole-assembly-design/candidate-whole-inputs.tsv'
ASSIGNMENTS = ROOT / 'results/ani-cases/2026-10-03/task12-astra-cp102-approved-final-reconciliation/accession-assignments.tsv'
DB = Path('/scratch/opendata/genomics/Kraken2/Refseq/Standard/20260226')
KRAKEN = '${USER_HOME}/sw/miniforge3/envs/kraken/bin/kraken2'


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, fields):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def fasta(path):
    header, sequence = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence).upper()
                header, sequence = line[1:], []
            else:
                if header is None:
                    raise ValueError('Sequence before FASTA header')
                sequence.append(line)
    if header is not None:
        yield header, ''.join(sequence).upper()


def tiles(sequence, length=150):
    for start in range(0, len(sequence), length):
        yield start, min(start + length, len(sequence)), sequence[start:start + length]


def ambiguous(sequence):
    return len(sequence) - sum(sequence.count(base) for base in 'ACGT')


def prepare(out, manifest=MANIFEST, assignments=ASSIGNMENTS, expected_count=1070):
    dest = out / 'inputs'
    dest.mkdir()  # Fail safely instead of overwriting a prior preparation.
    inputs, allrows = read(manifest), read(assignments)
    approved = {r['accession']: r for r in allrows if r['in_frozen_study'] == '1' and r['assignment_status'] != 'excluded_or_restricted_analysis_scope'}
    assert len(inputs) == expected_count == len(approved)
    assert len({r['accession'] for r in inputs}) == len(inputs)
    assert set(approved) == {r['accession'] for r in inputs}
    shutil.copy2(manifest, dest / 'candidate-whole-inputs.tsv')
    shutil.copy2(assignments, dest / 'accession-assignments.tsv')
    verified, contigs, accounting = [], [], []
    with open(dest / 'native.fna', 'w') as native, open(dest / 'fragments150.fna', 'w') as fragments, gzip.open(dest / 'fragment-map.tsv.gz', 'wt', compresslevel=1) as mapping:
        mapping.write('sequence_id\taccession\tcontig_id\tstart_0based\tend_exclusive\tlength\ttail\tambiguous_bases\tshorter_than_k35\n')
        for row in inputs:
            acc, path = row['accession'], Path(row['whole_fasta'])
            digest = sha(path)
            if digest != row['expected_sha256']:
                raise ValueError('Whole-assembly SHA256 mismatch: ' + acc)
            totals = collections.Counter()
            for index, (header, sequence) in enumerate(fasta(path), 1):
                if not sequence:
                    raise ValueError('Empty contig: ' + acc)
                cid = f'{acc}|c{index}'
                native.write(f'>{cid}\n{sequence}\n')
                contigs.append(dict(sequence_id=cid, accession=acc, source_header=header, source_sequence_accession=header.split()[0], length=len(sequence), ambiguous_bases=ambiguous(sequence)))
                totals.update(contigs=1, bases=len(sequence), ambiguous_bases=ambiguous(sequence))
                for start, end, fragment in tiles(sequence):
                    fid = f'{cid}|{start}|{end}'
                    tail, amb, short = int(end - start < 150), ambiguous(fragment), int(end - start < 35)
                    fragments.write(f'>{fid}\n{fragment}\n')
                    mapping.write(f'{fid}\t{acc}\t{cid}\t{start}\t{end}\t{end-start}\t{tail}\t{amb}\t{short}\n')
                    totals.update(fragments=1, fragment_bases=end-start, tail_fragments=tail, tail_bases=(end-start)*tail, short_k35_fragments=short, short_k35_bases=(end-start)*short, ambiguous_fragments=int(amb > 0))
            assert totals['contigs'] > 0 and totals['bases'] == totals['fragment_bases']
            verified.append(dict(row, observed_sha256=digest, hash_match=1, file_bytes=path.stat().st_size))
            accounting.append(dict(accession=acc, **{key: totals[key] for key in ['contigs', 'bases', 'ambiguous_bases', 'fragments', 'fragment_bases', 'tail_fragments', 'tail_bases', 'short_k35_fragments', 'short_k35_bases', 'ambiguous_fragments']}))
    write(dest / 'verified-inputs.tsv', verified, list(verified[0]))
    write(dest / 'contig-map.tsv', contigs, list(contigs[0]))
    write(dest / 'sequence-accounting.tsv', accounting, list(accounting[0]))
    write(dest / 'source-hashes.tsv', [dict(path=str(p), sha256=sha(p)) for p in [manifest, assignments, Path(__file__), ROOT/'scripts/apply_reporting_categories.py']], ['path', 'sha256'])
    (dest / 'PREPARED.json').write_text(json.dumps(dict(assemblies=len(inputs), totals=dict(sum((collections.Counter({k: int(v) for k, v in r.items() if k != 'accession'}) for r in accounting), collections.Counter()))), indent=2) + '\n')


def pin_database(out, db=DB):
    dest = out / 'database'
    dest.mkdir()
    records = []
    for name in ['hash.k2d', 'opts.k2d', 'taxo.k2d', 'nodes.dmp', 'names.dmp', 'seqid2taxid.map', 'library_report.tsv', 'ktaxonomy.tsv', 'inspect.txt']:
        path = db / name
        records.append(dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path)))
        if name != 'hash.k2d':
            shutil.copy2(path, dest / name)
    write(dest / 'database-hashes.tsv', records, ['path', 'bytes', 'sha256'])


def overlap(out):
    """Membership evidence is sequence/accession evidence, never independent-isolate truth."""
    db = out / 'database'
    seqids = set()
    with open(db / 'seqid2taxid.map') as handle:
        for line in handle:
            for token in re.findall(r'(?:NZ_)?[A-Z]+_[0-9]+\.[0-9]+|(?:NZ_)?[A-Z]{2}[0-9]+\.[0-9]+', line.split('\t')[0]):
                seqids.add(token)
    library = (db / 'library_report.tsv').read_text()
    assemblyids = set(re.findall(r'GC[AF]_[0-9]+\.[0-9]+', library))
    # Authoritative GCA/GCF pairing from the retained public NCBI snapshots.
    metadata_paths = [ROOT / p for p in [
        'results/ani-cases/2026-09-23/task12-astra-evidence/inputs/ncbi-2026-09-23.jsonl',
        'results/cohort/2026-09-15/task01-astra-public-provenance/assembly-reports.jsonl',
        'results/cohort/2026-09-15/task01-astra-public-provenance/routine-pairs/assembly-reports.jsonl',
        'results/cohort/2026-09-15/task01-astra-public-provenance/file-gaps/assembly-reports.jsonl',
        'results/ncbi/2026-09-15/file-gap-resolution/package/ncbi_dataset/data/assembly_data_report.jsonl',
    ]]
    # Keep this bounded to local public assembly metadata: no external downloads.
    pairs, pair_sources = {}, []
    for path in metadata_paths:
        try:
            with open(path) as handle:
                for line in handle:
                    obj = json.loads(line)
                    acc = obj.get('accession', '')
                    paired = (obj.get('assemblyInfo', {}).get('pairedAssembly', {}).get('accession', '')
                              or obj.get('assembly_info', {}).get('paired_assembly', {}).get('accession', '')
                              or obj.get('pairedAccession', '') or obj.get('paired_accession', ''))
                    if acc and paired:
                        pairs.setdefault(acc, set()).add(paired)
                        pairs.setdefault(paired, set()).add(acc)
            pair_sources.append(dict(path=str(path), sha256=sha(path)))
        except (ValueError, AttributeError):
            continue
    write(db / 'pairing-source-hashes.tsv', pair_sources, ['path', 'sha256'])
    rows, ctr = [], collections.defaultdict(collections.Counter)
    for r in read(out / 'inputs/contig-map.tsv'):
        seq = r['source_sequence_accession']
        exact = seq in seqids
        # NZ_ is an explicitly reported accession-wrapper sensitivity, not exact membership.
        wrapper = ('NZ_' + seq) in seqids or (seq.startswith('NZ_') and seq[3:] in seqids)
        rows.append(dict(accession=r['accession'], sequence_id=r['sequence_id'], source_sequence_accession=seq, exact_sequence_accession=int(exact), nz_wrapper_sequence_accession=int(wrapper)))
        ctr[r['accession']].update(contigs=1, exact_contigs=int(exact), exact_bases=int(r['length'])*exact, nz_wrapper_contigs=int(wrapper), nz_wrapper_bases=int(r['length'])*wrapper)
    write(db / 'sequence-overlap.tsv', rows, list(rows[0]))
    assemblyrows = []
    for acc, counts in sorted(ctr.items()):
        paired = sorted(pairs.get(acc, set()))
        assemblyrows.append(dict(accession=acc, paired_accessions=';'.join(paired), exact_assembly_accession_in_library=int(acc in assemblyids), paired_assembly_accessions_in_library=';'.join(p for p in paired if p in assemblyids), assembly_accession_audit='library_report_only_not_complete_build_inventory', **counts, independence='not_established_by_accession_membership_or_absence'))
    write(db / 'assembly-overlap.tsv', assemblyrows, list(assemblyrows[0]))


def load_taxonomy(db):
    parents, ranks, names = {}, {}, {}
    with open(db / 'nodes.dmp') as handle:
        for line in handle:
            fields = [v.strip() for v in line.split('|')]
            parents[fields[0]], ranks[fields[0]] = fields[1:3]
    with open(db / 'names.dmp') as handle:
        for line in handle:
            fields = [v.strip() for v in line.split('|')]
            if fields[3] == 'scientific name':
                names[fields[0]] = fields[1]
    return parents, ranks, names


def taxonomy_call(taxid, taxonomy):
    parents, ranks, names = taxonomy
    if taxid == '0':
        return dict(raw_taxid=taxid, raw_name='unclassified', raw_rank='unclassified', species_taxid='', species_name='', genus_taxid='', genus_name='', lineage_taxids='', resolution='unclassified')
    if taxid not in parents:
        raise ValueError('Raw taxid absent from pinned taxonomy: ' + taxid)
    lineage, seen, species, genus = [], set(), '', ''
    current = taxid
    while current not in seen:
        seen.add(current)
        lineage.append(current)
        if ranks.get(current) == 'species' and not species:
            species = current
        if ranks.get(current) == 'genus' and not genus:
            genus = current
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        current = parent
    return dict(raw_taxid=taxid, raw_name=names.get(taxid, ''), raw_rank=ranks.get(taxid, ''), species_taxid=species, species_name=names.get(species, ''), genus_taxid=genus, genus_name=names.get(genus, ''), lineage_taxids=';'.join(lineage), resolution='species_or_below' if species else 'genus_only' if ranks.get(taxid) == 'genus' else 'below_genus_without_named_species' if genus else 'above_genus_or_unresolved')


def compare_profile(row, call, cw):
    target = row['project_species']
    raw = call['raw_name']
    species = call['species_name']
    if call['raw_taxid'] == '0':
        return 'unclassified', 'no_call'
    if not target or target.endswith(' sp.'):
        return 'approved_unnamed_or_caution', 'organism_named_species_unresolved_no_error_inferred'
    if not species:
        return ('compatible_coarse_genus' if call['genus_name'] == target.split()[0] else 'unresolved_coarse_or_other_genus'), 'insufficient_species_resolution'
    parsed = parse_label(raw)
    groups = cw[1].get(parsed['name'], set()) if parsed['subspecies'] else set()
    targets = {s for s, _ in groups} if groups else cw[0].get(species, set())
    if not targets:
        return 'unmapped_taxonomy', 'no_validated_framework_mapping'
    if target not in targets:
        return 'framework_biological_discordance', 'sequence_call_outside_reviewed_organism_species_groups'
    if groups and row['project_subspecies'] and (target, row['project_subspecies']) not in groups and all(s for t, s in groups if t == target):
        return 'framework_subspecies_discordance', 'explicit_call_different_reviewed_subgroup'
    if species != target:
        return 'compatible_framework_name_or_rank', 'taxon_mapped_into_reviewed_species'
    return 'compatible_retained_species', 'same_species_or_coarser_within_species'


def parse_kraken(line, contigs, mode):
    fields = line.rstrip('\n').split('\t')
    if len(fields) != 5 or fields[0] not in {'C', 'U'}:
        raise ValueError('Malformed Kraken output')
    seqid, taxid, length = fields[1], fields[2], int(fields[3])
    parts = seqid.split('|')
    cid = '|'.join(parts[:2])
    contig = contigs[cid]
    if mode == 'native':
        if len(parts) != 2 or length != int(contig['length']):
            raise ValueError('Native contig length/identity mismatch')
    else:
        if len(parts) != 4:
            raise ValueError('Fragment ID lacks coordinates')
        start, end = map(int, parts[2:])
        if start % 150 or end != min(start + 150, int(contig['length'])) or length != end - start or not 0 <= start < end <= int(contig['length']):
            raise ValueError('Fragment coordinates/length invalid')
    if (fields[0] == 'U') != (taxid == '0'):
        raise ValueError('Inconsistent classification status and taxid')
    return contig['accession'], cid, taxid, length


def summarize(out):
    contigs = {r['sequence_id']: r for r in read(out / 'inputs/contig-map.tsv')}
    accounting = {r['accession']: r for r in read(out / 'inputs/sequence-accounting.tsv')}
    assignments = read(out / 'inputs/accession-assignments.tsv')
    byacc, cw = {r['accession']: r for r in assignments}, make_crosswalk(assignments)
    taxonomy, calls = load_taxonomy(out / 'database'), {}
    profiles, genome_summaries, category_profiles, validations = [], [], [], []
    for mode in ['native', 'fragments150']:
        for confidence in ['0', '0.1']:
            counters, seen_contigs = collections.defaultdict(collections.Counter), collections.Counter()
            # Sum and sum-of-squares of tile indices catch missing/duplicate IDs cheaply;
            # output record count and lengths also checked against every source contig.
            idsums = collections.defaultdict(collections.Counter)
            seen_tiles = {cid: bytearray((1 if mode == 'native' else (int(r['length'])+149)//150)) for cid, r in contigs.items()}
            path = out / 'raw' / f'{mode}.confidence{confidence}.output.gz'
            with gzip.open(path, 'rt') as handle:
                for line in handle:
                    acc, cid, taxid, length = parse_kraken(line, contigs, mode)
                    counters[(acc, taxid)].update(sequences=1, bases=length)
                    seen_contigs[cid] += 1
                    index = 0 if mode == 'native' else int(line.split('\t')[1].split('|')[2]) // 150
                    if seen_tiles[cid][index]:
                        raise ValueError('Duplicate Kraken sequence ID: '+line.split('\t')[1])
                    seen_tiles[cid][index] = 1
                    idsums[cid].update(index_sum=index, index_square_sum=index*index, bases=length)
                    if taxid not in calls:
                        calls[taxid] = taxonomy_call(taxid, taxonomy)
            for cid, contig in contigs.items():
                n = 1 if mode == 'native' else (int(contig['length']) + 149)//150
                assert seen_contigs[cid] == n, (cid, seen_contigs[cid], n)
                assert idsums[cid]['bases'] == int(contig['length'])
                assert idsums[cid]['index_sum'] == n*(n-1)//2
                assert idsums[cid]['index_square_sum'] == n*(n-1)*(2*n-1)//6
            totals = collections.defaultdict(collections.Counter)
            categories = collections.defaultdict(collections.Counter)
            for (acc, taxid), count in sorted(counters.items()):
                a, row, call = accounting[acc], byacc[acc], calls[taxid]
                denom = int(a['contigs'] if mode == 'native' else a['fragments'])
                category, reason = compare_profile(row, call, cw)
                profiles.append(dict(accession=acc, input_mode=mode, confidence=confidence, **call, **count, sequence_denominator=denom, base_denominator=a['bases'], sequence_fraction=count['sequences']/denom, represented_base_fraction=count['bases']/int(a['bases']), project_species=row['project_species'], project_subspecies=row['project_subspecies'], approved_flags=row['nomenclature_and_boundary_flags'], literal_species_name_difference=int(bool(call['species_name'] and call['species_name'] != row['project_species'])), framework_comparison=category, comparison_reason=reason))
                totals[acc].update(count)
                categories[(acc, category)].update(count)
            for acc, total in sorted(totals.items()):
                a = accounting[acc]
                assert total['bases'] == int(a['bases'])
                expected = int(a['contigs'] if mode == 'native' else a['fragments'])
                assert total['sequences'] == expected
                unclassified = counters.get((acc, '0'), collections.Counter())
                genome_summaries.append(dict(accession=acc, input_mode=mode, confidence=confidence, sequences=total['sequences'], bases=total['bases'], unclassified_sequences=unclassified['sequences'], unclassified_bases=unclassified['bases'], unclassified_sequence_fraction=unclassified['sequences']/expected, unclassified_base_fraction=unclassified['bases']/total['bases'], hard_genome_call='not_defined_profile_only', **{k: v for k, v in a.items() if k not in {'accession', 'bases'}}))
            for (acc, category), count in sorted(categories.items()):
                denom = int(accounting[acc]['contigs'] if mode == 'native' else accounting[acc]['fragments'])
                category_profiles.append(dict(accession=acc, input_mode=mode, confidence=confidence, category=category, **count, sequence_fraction=count['sequences']/denom, represented_base_fraction=count['bases']/int(accounting[acc]['bases'])))
            assert set(totals) == set(accounting)
            validations.append(dict(mode=mode, confidence=confidence, assemblies=len(totals), sequences=sum(v['sequences'] for v in totals.values()), bases=sum(v['bases'] for v in totals.values()), status='PASS'))
    dest = out / 'summaries'
    dest.mkdir()
    for name, rows in [('taxid-profiles.tsv', profiles), ('assembly-accounting.tsv', genome_summaries), ('framework-category-profiles.tsv', category_profiles), ('taxonomy-calls.tsv', list(calls.values()))]:
        write(dest / name, rows, list(rows[0]))
    (dest / 'validation.json').write_text(json.dumps(validations, indent=2) + '\n')


def run_command(command, out, name):
    start = time.time()
    with open(out / 'logs' / (name + '.stdout'), 'w') as stdout, open(out / 'logs' / (name + '.stderr'), 'w') as stderr:
        result = subprocess.run(['/usr/bin/time', '-v', '-o', str(out / 'logs' / (name + '.time')), *command], stdout=stdout, stderr=stderr)
    with open(out / 'commands.jsonl', 'a') as handle:
        handle.write(json.dumps(dict(command=command, name=name, start_unix=start, end_unix=time.time(), exit_status=result.returncode)) + '\n')
    if result.returncode:
        raise RuntimeError(f'{name} failed with {result.returncode}')


def pipeline(out, threads):
    start = time.time()
    version = subprocess.check_output([KRAKEN, '--version'], text=True)
    if 'version 2.1.6' not in version:
        raise ValueError('Unexpected Kraken version')
    (out / 'kraken-version.txt').write_text(version)
    for name, action in [('database-pin', lambda: pin_database(out)), ('prepare', lambda: prepare(out)), ('overlap', lambda: overlap(out))]:
        began = time.time()
        action()
        with open(out / 'stages.jsonl', 'a') as handle:
            handle.write(json.dumps(dict(stage=name, elapsed_seconds=time.time()-began, exit_status=0)) + '\n')
    raw = out / 'raw'
    raw.mkdir()
    for mode in ['native', 'fragments150']:
        for confidence in ['0', '0.1']:
            name = f'{mode}.confidence{confidence}'
            output = raw / (name + '.output')
            command = [KRAKEN, '--db', str(DB), '--threads', str(threads), '--confidence', confidence, '--report', str(raw / (name + '.report.tsv')), '--output', str(output), str(out / 'inputs' / (mode + '.fna'))]
            run_command(command, out, name)
            run_command(['gzip', '-1', str(output)], out, name + '.compress')
    began = time.time()
    summarize(out)
    with open(out / 'stages.jsonl', 'a') as handle:
        handle.write(json.dumps(dict(stage='summarize', elapsed_seconds=time.time()-began, exit_status=0)) + '\n')
    hashes = []
    for path in sorted(out.rglob('*')):
        if path.is_file() and path.name != 'output-hashes.tsv':
            hashes.append(dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size))
    write(out / 'output-hashes.tsv', hashes, ['path', 'sha256', 'bytes'])
    (out / 'COMPLETED.json').write_text(json.dumps(dict(exit_status=0, elapsed_seconds=time.time()-start, summary_validation='summaries/validation.json'), indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['pipeline', 'prepare', 'summarize', 'overlap'])
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'logs').mkdir(exist_ok=True)
    if args.operation == 'pipeline':
        pipeline(args.output.resolve(), args.threads)
    else:
        globals()[args.operation](args.output.resolve())


if __name__ == '__main__':
    main()
