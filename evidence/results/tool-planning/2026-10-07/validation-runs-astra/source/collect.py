"""Account for every attempted genome and verify immutable native outputs before scoring."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from common import RUN, ROOT, read, write, sha, json_write, verify_object, command
csv.field_size_limit(1_000_000_000)


def audit_query(row, root=RUN):
    directory = root / 'queries' / row['coded_id'] / 'attempt01'
    record = json.loads((directory / 'finished.json').read_text())  # Missing is not a failed biological run.
    if (record['genome_id'] != row['genome_id'] or record['coded_id'] != row['coded_id']
        or record['query_sha256'] != row['query_sha256'] or record['query_unchanged'] is not True
        or record['infrastructure_error'] or not record['process_started']
        or type(record['exit_status']) is not int):
        raise ValueError('unresolved execution/input failure: ' + row['coded_id'])
    if record['argv'] != command(Path(row['fasta_path']), directory / 'analyze'):
        raise ValueError('native command differs from frozen command')
    item = dict(genome_id=row['genome_id'], coded_id=row['coded_id'], layer=row['layer'],
                exit_status=record['exit_status'], seconds=record['seconds'], near_clone='not_assessed',
                child_peak_rss_kib='', report_tsv='', report_sha256='', failure_reason='',
                finished_sha256=sha(directory / 'finished.json'), native_run_sha256='')
    if record['exit_status'] != 0:
        failure = directory / 'analyze/failure.json'
        item['failure_reason'] = json.loads(failure.read_text()).get('error', '') if failure.exists() else 'see native stderr'
        return item
    output = directory / 'analyze'
    run = json.loads((output / 'run.json').read_text())
    if run['status'] != 'complete' or run['positive_calling_policy_enabled'] is not False:
        raise ValueError('unexpected native completion/policy status')
    inventory = {str(p.relative_to(output)) for p in output.rglob('*') if p.is_file() and p != output / 'run.json'}
    if set(run['artifacts_sha256']) != inventory or not {'evidence.json', 'report.tsv'} <= inventory:
        raise ValueError('native artifact manifest incomplete')
    for name, digest in run['artifacts_sha256'].items():
        artifact = (output / name).resolve()
        if not artifact.is_relative_to(output.resolve()) or sha(artifact) != digest:
            raise ValueError('native artifact changed: ' + row['coded_id'] + '/' + name)
    report_path = output / 'report.tsv'
    if record['report_tsv'] != str(report_path) or sha(report_path) != record['report_sha256']:
        raise ValueError('finished report path/digest mismatch')
    reports = read(report_path)
    if len(reports) != 1:
        raise ValueError('one genome must produce one report row')
    report = reports[0]
    if report['query_sha256'] != row['query_sha256'] or report['assigned_species_id'] or report['assigned_subspecies_id']:
        raise ValueError('query/positive assignment violation')
    evidence = json.loads((output / 'evidence.json').read_text())
    near_clone = any(h['ani_percent'] >= 99.9 and h['aligned_fraction_query_percent'] >= 90
                     and h['aligned_fraction_reference_percent'] >= 90 for h in evidence['ranked_reference_evidence'])
    item.update(near_clone=str(near_clone), child_peak_rss_kib=run['child_peak_rss_kib'],
                report_tsv=str(report_path), report_sha256=sha(report_path), native_run_sha256=sha(output / 'run.json'))
    return item


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--layer', choices=['layer1', 'layer3'], required=True)
    args = p.parse_args()
    verify_object()
    for pin in read(RUN / 'inputs/execution-pins.tsv') + read(RUN / 'inputs/driver-pins.tsv'):
        if sha(pin['path']) != pin['sha256']:
            raise ValueError('execution input/driver changed')
    rows = [r for r in read(RUN / 'inputs/execution-manifest.tsv') if r['layer'] == args.layer]
    if len({r['genome_id'] for r in rows}) != len(rows):
        raise ValueError('duplicate execution identity')
    for index in range((len(rows)+19)//20):
        worker = json.loads((RUN / 'logs' / f'{args.layer}-worker-{index:04d}/finished.json').read_text())
        if worker['infrastructure_error'] or worker['incomplete_coded_ids'] or worker['unattempted_coded_ids']:
            raise ValueError('incomplete worker: ' + str(index))
    ledger = [audit_query(row) for row in rows]
    output = RUN / 'outputs' / (args.layer + '-collection')
    output.mkdir(exist_ok=False)
    write(output / 'attempt-ledger.tsv', ledger)
    successful = [r for r in ledger if r['exit_status'] == 0]
    write(output / 'results.tsv', [{k: r[k] for k in ('genome_id', 'report_tsv', 'report_sha256')}
                                   for r in successful], ['genome_id', 'report_tsv', 'report_sha256'])
    summary = dict(layer=args.layer, genomes=len(rows), successful=len(successful),
                   failed=len(rows)-len(successful), all_accounted=True,
                   near_clones=dict(Counter(r['near_clone'] for r in ledger)),
                   native_seconds=sum(r['seconds'] for r in ledger), new_analysis=False,
                   production_assignments_enabled=False, collector_sha256=sha(__file__))
    json_write(output / 'summary.json', summary)
    print(json.dumps(summary))

if __name__ == '__main__':
    main()
