#!/usr/bin/env python3
"""Independently audit raw reciprocal pyANI BLAST outputs and scope inputs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def raw_summary(path):
    """Apply documented pyANI filters, first passing HSP per query fragment."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing raw BLAST evidence: {path}')
    chosen = {}
    rows = 0
    with path.open() as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if len(row) != 15:
                raise ValueError('Expected the pinned 15-column BLAST format')
            rows += 1
            q, _, length, mismatch, pid, nident, qlen, *_ = row
            length, mismatch, qlen, gaps = map(int, (length, mismatch, qlen, row[14]))
            ungapped = length - gaps
            if ungapped / qlen > .7 and (ungapped - mismatch) / qlen > .3:
                chosen.setdefault(q, (float(pid), ungapped, mismatch + gaps))
    if not chosen:
        raise ValueError('No passing raw BLAST hits')
    return dict(raw_hsps=rows, passing_fragments=len(chosen),
                ani_percent=sum(x[0] for x in chosen.values()) / len(chosen),
                aligned_query_bases=sum(x[1] for x in chosen.values()),
                similarity_errors=sum(x[2] for x in chosen.values()))


def matrix(path):
    with path.open() as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    return {(r[0], c): float(v) for r in rows[1:] for c, v in zip(rows[0][1:], r[1:])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    run = args.run
    records = json.loads((run / 'inputs/manifest.json').read_text())
    assert len(records) == 4
    for row in records:
        assert hashlib.sha256(Path(row['path']).read_bytes()).hexdigest() == row['sha256']
    with (run / 'inputs/stats.tsv').open() as handle:
        stats = list(csv.DictReader(handle, delimiter='\t'))
    lens = {(Path(r['file']).parent.name, Path(r['file']).stem): int(r['sum_len']) for r in stats}
    output = []
    for scope in ['chromosome', 'all_replicons']:
        out = run / ('anib-' + scope + '-v2')
        cmd = json.loads((run / 'logs' / (scope + '-v2.command.json')).read_text())
        assert cmd['exit_code'] == 0
        assert 'At least one BLAST run failed' not in (run / 'logs' / (scope + '-v2.pyani.log')).read_text()
        matrices = {n: matrix(out / ('ANIb_' + n + '.tab')) for n in
                    ['percentage_identity', 'alignment_coverage', 'alignment_lengths']}
        for q, s in [('GCA_000025565.1', 'GCA_009660455.1'), ('GCA_009660455.1', 'GCA_000025565.1')]:
            raw = out / 'blastn_output' / (q + '_vs_' + s + '.blast_tab')
            result = raw_summary(raw)
            coverage = result['aligned_query_bases'] / lens[scope, q] * 100
            for name, value in [('percentage_identity', result['ani_percent'] / 100),
                                ('alignment_coverage', coverage / 100),
                                ('alignment_lengths', result['aligned_query_bases'])]:
                assert math.isclose(matrices[name][q, s], value, abs_tol=1e-10), (scope, q, name)
            assert 0 < coverage <= 100 and 0 < result['ani_percent'] < 100
            output.append(dict(scope=scope, query=q, subject=s, **result,
                               query_length=lens[scope, q], query_coverage_percent=coverage,
                               raw_path=str(raw), raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest()))
    with (run / 'directional-results.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(output)
    report = dict(status='PASS', input_hashes=4, raw_directions=4, matrix_values_recomputed=12,
                  similarity_errors_matrix_not_used='Installed pyANI symmetrizes/overwrites directional error counts; raw totals retained',
                  results=output)
    scoring = []
    for scope in ['chromosome', 'all_replicons']:
        folder = run / ('scoring-' + scope)
        commands = json.loads((folder / 'commands.json').read_text())
        assert len(commands) == 2 and all(c['exit_code'] == 0 for c in commands)
        for q, s in [('GCA_000025565.1', 'GCA_009660455.1'), ('GCA_009660455.1', 'GCA_000025565.1')]:
            raw = folder / (q + '_vs_' + s + '.blast_tab')
            result = raw_summary(raw)
            scoring.append(dict(scope=scope, query=q, subject=s, **result,
                                query_length=lens[scope,q],
                                query_coverage_percent=100*result['aligned_query_bases']/lens[scope,q]))
    with (run / 'scoring-sensitivity.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scoring[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(scoring)
    sets = {}
    with (run / 'inputs/record-hashes.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (Path(row['source_file']).parent.name, Path(row['source_file']).stem)
            sets.setdefault(key, {})[row['#id']] = (int(row['length']), row['seq.hash'])
    assert sets['chromosome', 'GCA_009660455.1'] == sets['all_replicons', 'GCA_009660455.1']
    assert set(sets['chromosome', 'GCA_009660455.1']) == {f'WJWQ01{i:06d}.1' for i in range(1,52)}
    assert sets['chromosome', 'GCA_000025565.1']['CP001918.1'] == sets['all_replicons', 'GCA_000025565.1']['CP001918.1']
    assert set(sets['all_replicons', 'GCA_000025565.1']) == {'CP001918.1', 'CP001919.1', 'CP001920.1'}
    genbank = (run / 'sources/CP001918.1.gb').read_text().split('ORIGIN', 1)[1].split('//', 1)[0]
    sequence = ''.join(c for c in genbank if c.isalpha()).lower()
    assert len(sequence) == 5314581
    assert hashlib.md5(sequence.encode()).hexdigest() == sets['chromosome', 'GCA_000025565.1']['CP001918.1'][1]
    report.update(scoring_directions=4, scope_checks=dict(matched_dissolvens_records=51,
                  matched_cloacae_chromosome_records=1, public_CP001918_1_sequence_matches=True,
                  dissolvens_wgs_accessions_complete=True, all_replicon_accessions_explicit=True))
    (run / 'validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
