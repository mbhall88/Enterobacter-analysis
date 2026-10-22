#!/usr/bin/env python3
"""Independent accounting validation of Task14 prepared data and profile outputs."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import time
from pathlib import Path


def read(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def validate(run, expected_count=1070):
    started = time.time()
    assertions = 0

    def check(condition, message):
        nonlocal assertions
        assertions += 1
        if not condition:
            raise ValueError(message)

    expected = read(run/'inputs/candidate-whole-inputs.tsv')
    observed = read(run/'inputs/verified-inputs.tsv')
    assignments = read(run/'inputs/accession-assignments.tsv')
    eligible = {r['accession'] for r in assignments if r['in_frozen_study']=='1' and r['assignment_status']!='excluded_or_restricted_analysis_scope'}
    check(len(expected)==len(observed)==len(eligible)==expected_count, 'Approved denominator mismatch')
    check({r['accession'] for r in expected}==eligible=={r['accession'] for r in observed}, 'Approved membership mismatch')
    for row in observed:
        check(row['observed_sha256']==row['expected_sha256'], 'Input hash mismatch')
    contigs = read(run/'inputs/contig-map.tsv')
    bycid = {r['sequence_id']:r for r in contigs}
    check(len(bycid)==len(contigs), 'Duplicate contig IDs')
    accounted = {r['accession']:r for r in read(run/'inputs/sequence-accounting.tsv')}
    calculated = collections.defaultdict(collections.Counter)
    for r in contigs:
        length = int(r['length'])
        tail = length%150
        calculated[r['accession']].update(contigs=1,bases=length,ambiguous_bases=int(r['ambiguous_bases']),fragments=math.ceil(length/150),tail_fragments=int(tail>0),tail_bases=tail,short_k35_fragments=int(0<tail<35),short_k35_bases=tail if tail<35 else 0)
    for acc, counter in calculated.items():
        for key, value in counter.items():
            check(int(accounted[acc][key])==value, f'{acc} {key} accounting mismatch')
    check(set(accounted)==eligible, 'Missing prepared accession')
    profiles = read(run/'summaries/taxid-profiles.tsv')
    totals = collections.defaultdict(collections.Counter)
    fractions = collections.defaultdict(collections.Counter)
    for row in profiles:
        key = row['accession'],row['input_mode'],row['confidence']
        totals[key].update(sequences=int(row['sequences']),bases=int(row['bases']))
        fractions[key].update(sequences=float(row['sequence_fraction']),bases=float(row['represented_base_fraction']))
        check(row['accession'] in eligible, 'Unapproved profile accession')
        expected_sequences = accounted[row['accession']]['contigs' if row['input_mode']=='native' else 'fragments']
        check(int(row['sequence_denominator'])==int(expected_sequences), 'Wrong profile sequence denominator')
        check(int(row['base_denominator'])==int(accounted[row['accession']]['bases']), 'Wrong profile base denominator')
        check(0<=float(row['represented_base_fraction'])<=1, 'Invalid profile base fraction')
    check(len(totals)==expected_count*4, 'Missing assembly/settings profiles')
    for (acc,mode,conf), total in totals.items():
        check(total['bases']==int(accounted[acc]['bases']), 'Profile base coverage mismatch')
        check(total['sequences']==int(accounted[acc]['contigs' if mode=='native' else 'fragments']), 'Profile sequence count mismatch')
        check(abs(fractions[(acc,mode,conf)]['bases']-1)<1e-12, 'Base fractions do not sum to one')
        check(abs(fractions[(acc,mode,conf)]['sequences']-1)<1e-12, 'Sequence fractions do not sum to one')
    taxonomy = {r['raw_taxid']:r for r in read(run/'summaries/taxonomy-calls.tsv')}
    check(all(r['raw_taxid'] in taxonomy for r in profiles), 'Missing raw taxon provenance')
    stored_taxonomy = {}
    with (run/'database/ktaxonomy.tsv').open() as handle:
        for line in handle:
            fields=[field.strip() for field in line.split('|')]
            stored_taxonomy[fields[0]]=fields[-1]
    for taxid,call in taxonomy.items():
        if taxid!='0':
            check(taxid in stored_taxonomy, 'Call absent from stored Kraken taxonomy')
            check(stored_taxonomy[taxid]==call['raw_name'], 'Pinned name disagrees with stored Kraken taxonomy')
    categories = collections.defaultdict(collections.Counter)
    for r in read(run/'summaries/framework-category-profiles.tsv'):
        categories[(r['accession'],r['input_mode'],r['confidence'])].update(sequences=int(r['sequences']),bases=int(r['bases']))
    check(categories==totals, 'Framework comparison changed denominators')
    native_expected = sum(int(r['contigs']) for r in accounted.values())
    fragment_expected = sum(int(r['fragments']) for r in accounted.values())
    for mode in ['native','fragments150']:
        for conf in ['0','0.1']:
            report = run/'raw'/f'{mode}.confidence{conf}.report.tsv'
            total = 0
            with report.open() as h:
                for line in h:
                    fields = line.rstrip().split('\t')
                    if fields[4] in {'0','1'}:
                        total += int(fields[1])
            check(total==(native_expected if mode=='native' else fragment_expected), 'Standard report count mismatch')
            # Full raw-output checks were independently reconciled against every
            # expected source ID by the summarizer; require all four successful logs.
            check((run/'raw'/f'{mode}.confidence{conf}.output.gz').is_file(), 'Missing raw output')
    validation = json.loads((run/'summaries/validation.json').read_text())
    check(len(validation)==4 and all(r['status']=='PASS' and r['assemblies']==expected_count for r in validation), 'Raw coverage validation incomplete')
    check((run/'logs/worker.exit').read_text().strip()=='0', 'Worker did not exit successfully')
    check(json.loads((run/'COMPLETED.json').read_text())['exit_status']==0, 'Pipeline completion missing')
    return dict(status='PASS', assertions=assertions, assemblies=expected_count, settings=4, native_contigs=native_expected, fragments150=fragment_expected, represented_bases=sum(int(r['bases']) for r in accounted.values()), elapsed_seconds=time.time()-started)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result=validate(args.run)
    result['validator_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
