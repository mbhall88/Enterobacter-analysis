#!/usr/bin/env python3
"""Audit both saved NUCmer directions using the pinned pyANI parser; no new alignment.

Run with .work/envs/task12-pyani/bin/python. The legacy matrix stores the two
identities, but its lengths/coverage cells are overwritten by the last sorted
direction. Retain all four coverages explicitly and verify the exported matrix.
"""
import argparse
import json
from pathlib import Path
from pyani.anim import parse_delta
from audit_ani_tables import sha256_file
from build_case_register import read_tsv, write_tsv, passes
from summarize_task12_sensitivities import matrix


def audit(run, output):
    output.mkdir(parents=True, exist_ok=False)
    results = []
    comparisons = []
    for pair in read_tsv(run / 'anim/selected-pairs.tsv'):
        a, b = pair['a'], pair['b']
        work = run / 'anim/pairs' / (a + '__' + b) / 'work'
        lengths = matrix(work / 'ANIm_alignment_lengths.tab')
        cover = matrix(work / 'ANIm_alignment_coverage.tab')
        identity = matrix(work / 'ANIm_percentage_identity.tab')
        filters = sorted((work / 'nucmer_output').glob('*/*.filter'))
        if len(filters) != 2:
            raise ValueError(f'Expected both NUCmer directions: {a}, {b}')
        observed = []
        for path in filters:
            ref, query = path.stem.split('_vs_')
            lr, lq, ani, errors = parse_delta(path)
            if abs(identity[ref][query] - ani) > 1e-12:
                raise ValueError('Identity matrix and raw alignment disagree')
            row = dict(a=a, b=b, nucmer_reference=ref, nucmer_query=query,
                       identity=100*ani, reference_aligned_bases=lr,
                       query_aligned_bases=lq, reference_total_bases=lengths[ref][ref],
                       query_total_bases=lengths[query][query],
                       reference_coverage=100*lr/lengths[ref][ref],
                       query_coverage=100*lq/lengths[query][query],
                       similarity_errors=errors, source=str(path), sha256=sha256_file(path))
            results.append(row)
            observed.append(row)
        # pyANI overwrites both length/coverage cells in sorted delta-file order.
        last = observed[-1]
        ref, query = last['nucmer_reference'], last['nucmer_query']
        for acc, field in [(ref, 'reference'), (query, 'query')]:
            other = query if acc == ref else ref
            if abs(cover[acc][other]*100-last[field+'_coverage']) > 1e-10:
                raise ValueError('Coverage matrix does not match its documented last-direction semantics')
            if lengths[acc][other] != last[field+'_aligned_bases']:
                raise ValueError('Length matrix and raw alignment disagree')
        min_ani = min(r['identity'] for r in observed)
        min_coverage = min(r[k] for r in observed for k in ['reference_coverage', 'query_coverage'])
        matrix_rule = passes(min_ani, cover[a][b]*100, cover[b][a]*100)
        all_direction_rule = min_ani >= 95 and min_coverage >= 80
        comparisons.append(dict(a=a,b=b,minimum_of_two_identities=min_ani,
                                minimum_of_four_coverages=min_coverage,
                                matrix_rule_passes=matrix_rule,
                                all_directions_rule_passes=all_direction_rule,
                                matrix_coverage_source=str(filters[-1])))
    write_tsv(output / 'directional-alignments.tsv', results)
    write_tsv(output / 'matrix-semantics-audit.tsv', comparisons)
    (output / 'summary.json').write_text(json.dumps(dict(
        pairs=len(comparisons), directional_alignments=len(results),
        matrix_checks_passed=True,
        all_direction_rule_outcome_changes=sum(r['matrix_rule_passes'] != r['all_directions_rule_passes'] for r in comparisons),
        coverage_semantics='Both matrix coverage cells come from the last sorted NUCmer direction; all four directional coverages retained separately.'
    ), indent=2)+'\n')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    audit(a.run,a.output)
