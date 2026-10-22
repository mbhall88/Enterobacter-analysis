#!/usr/bin/env python3
"""Derive versioned informal 99.5%-ANI groups from the frozen Group A ledger.

No ANI is measured here. Complete linkage is a project implementation choice;
the cited paper supplies the threshold, not a mandated linkage algorithm.
"""
import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, fields=None):
    rows = list(rows)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def partition(ani, ids, method='complete', cutoff=99.5):
    """Sort accessions before linkage and number by minimum accession in group."""
    ids = np.asarray(ids)
    order = np.argsort(ids)
    if len(set(ids)) != len(ids) or ani.shape != (len(ids), len(ids)):
        raise ValueError('Nonunique IDs or invalid shape')
    if not (np.isfinite(ani).all() and np.allclose(ani, ani.T)
            and np.all(np.diag(ani) == 100) and np.all((ani >= 0) & (ani <= 100))):
        raise ValueError('Invalid ANI matrix')
    if len(ids) == 1:
        return np.array(['GV01'])
    z = linkage(squareform(100 - ani[np.ix_(order, order)]), method=method)
    raw = fcluster(z, 100 - cutoff + 1e-10, criterion='distance')
    groups = sorted((sorted(ids[order][raw == k]) for k in set(raw)), key=lambda g: g[0])
    mapping = {accession: f'GV{i:02d}' for i, g in enumerate(groups, 1) for accession in g}
    return np.array([mapping[a] for a in ids])


def counts(values):
    return ';'.join(f'{k}:{v}' for k, v in sorted(Counter(values).items()))


def run(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    assignments = read(args.populations / 'outputs/assignments.tsv')
    ids = sorted(r['accession'] for r in assignments)
    n = len(ids)
    ix = {a: i for i, a in enumerate(ids)}
    ani = np.full((n, n), np.nan)
    np.fill_diagonal(ani, 100)
    pairs = read(args.populations / 'inputs/pairs.tsv')
    checks = []

    def check(name, condition):
        checks.append({'check': name, 'pass': bool(condition)})
        if not condition:
            raise AssertionError(name)

    check('173 unique genomes and all 14878 unordered pairs', n == 173 and len(ix) == n and len(pairs) == n*(n-1)//2)
    for r in pairs:
        i, j = ix[r['a']], ix[r['b']]
        check(f'unique measured chromosome skani pair {i}-{j}',
              i != j and np.isnan(ani[i, j]) and r['status'] == 'measured'
              and r['method'] == 'skani_0.3.1' and r['scope'] == 'chromosome_to_chromosome'
              and min(float(r['af_a']), float(r['af_b'])) >= 80)
        ani[i, j] = ani[j, i] = float(r['ani'])
    # Independently reconcile TSV against the previous compressed matrix.
    old = np.load(args.populations / 'inputs/matrices.npz')
    oldix = {str(a): i for i, a in enumerate(old['ids'])}
    order = [oldix[a] for a in ids]
    check('TSV matches frozen matrix', np.array_equal(ani, old['ani'][np.ix_(order, order)]))
    overlay = {r['accession']: r for r in read(args.naming)}
    check('overlay exact membership', set(overlay) == set(ids))
    labels = partition(ani, ids)
    methods = {m: partition(ani, ids, m) for m in ['complete', 'average', 'single']}
    membership, summaries, sensitivity, threshold_pairs = [], [], [], []
    for i, a in enumerate(ids):
        row = dict(overlay[a])
        row.update(genomovar_id=labels[i], genomovar_status='informal_operational_ANI_group',
                   genomovar_rule='skani_chromosome_ANI_ge_99.5_complete_linkage',
                   genomovar_version='2026-10-02-v1',
                   average_linkage_group=methods['average'][i], single_linkage_group=methods['single'][i])
        membership.append(row)
    for gv in sorted(set(labels)):
        members = np.where(labels == gv)[0]
        values = ani[np.ix_(members, members)][np.triu_indices(len(members), 1)]
        check(f'{gv} all internal pairs meet threshold', np.all(values >= 99.5 - 1e-10))
        # ANI medoid: maximum mean ANI to the other members; accession breaks ties.
        means = ani[np.ix_(members, members)].mean(axis=1)
        representative = ids[members[int(np.argmax(means))]]
        rows = [overlay[ids[i]] for i in members]
        summaries.append(dict(genomovar_id=gv, genomes=len(members),
                              representative_accession=representative,
                              within_ANI_min='' if len(values) == 0 else f'{values.min():.2f}',
                              within_ANI_max='' if len(values) == 0 else f'{values.max():.2f}',
                              population_counts=counts(r['operational_population'] for r in rows),
                              exact_ST_counts=counts(r['ST'] for r in rows if r['ST_status'] == 'exact_ST'),
                              ST_unassigned=sum(r['ST_status'] != 'exact_ST' for r in rows),
                              historical_types=';'.join(r['historical_type_relationship'] for r in rows if r['historical_type_relationship'])))
    for method, cutoff in itertools.product(methods, [99.495, 99.5, 99.505]):
        ll = partition(ani, ids, method, cutoff)
        cc = Counter(ll)
        coassign = ll[:, None] == ll[None, :]
        sensitivity.append(dict(linkage=method, cutoff=cutoff, groups=len(cc),
                                singleton_groups=sum(v == 1 for v in cc.values()),
                                largest_group=max(cc.values()),
                                same_partition_as_primary=bool(np.array_equal(coassign, labels[:, None] == labels[None, :])),
                                within_pairs_below_99_5=int(np.sum(np.triu(coassign & (ani < 99.5), 1)))))
    for i, j in itertools.combinations(range(n), 2):
        if ani[i, j] >= 99.5 and labels[i] != labels[j] or abs(ani[i, j] - 99.5) < 1e-8:
            threshold_pairs.append(dict(a=ids[i], b=ids[j], ANI=ani[i, j], a_group=labels[i], b_group=labels[j],
                                        exactly_at_cutoff=bool(abs(ani[i, j] - 99.5) < 1e-8)))
    psummary = []
    for pop in sorted({r['operational_population'] for r in membership}):
        rr = [r for r in membership if r['operational_population'] == pop]
        psummary.append(dict(population=pop, genomes=len(rr), genomovars=len({r['genomovar_id'] for r in rr}),
                             membership=counts(r['genomovar_id'] for r in rr)))
    check('all genomes assigned exactly once', len(membership) == n and sum(s['genomes'] for s in summaries) == n)
    check('ST calls preserved', sum(r['ST_status'] == 'exact_ST' for r in membership) == 168)
    check('P3 intact', len({r['genomovar_id'] for r in membership if r['operational_population'] == 'P3'}) == 1)
    write(out / 'memberships.tsv', membership)
    write(out / 'genomovars.tsv', summaries)
    write(out / 'population-crosswalk.tsv', psummary)
    write(out / 'sensitivity.tsv', sensitivity)
    write(out / 'threshold-bridge-pairs.tsv', threshold_pairs,
          ['a', 'b', 'ANI', 'a_group', 'b_group', 'exactly_at_cutoff'])
    write(out / 'validation.tsv', checks)
    summary = dict(genomes=n, pairs=len(pairs), groups=len(summaries),
                   singleton_groups=sum(s['genomes'] == 1 for s in summaries),
                   group_sizes=sorted([s['genomes'] for s in summaries], reverse=True),
                   checks_passed=len(checks), sensitivity=sensitivity,
                   versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    inputs = [args.populations / f for f in ['outputs/assignments.tsv', 'inputs/pairs.tsv', 'inputs/matrices.npz']] + [args.naming, Path(__file__)]
    (out.parent / 'analysis-inputs.sha256').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p}\n' for p in inputs))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--populations', type=Path, required=True)
    parser.add_argument('--naming', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
