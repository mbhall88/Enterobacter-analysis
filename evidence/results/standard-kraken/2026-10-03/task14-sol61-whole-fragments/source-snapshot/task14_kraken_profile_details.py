#!/usr/bin/env python3
"""Compact native-contig calls and fragment-quality strata from preserved raw output."""
import argparse
import collections
import gzip
import json
import statistics
import time
from pathlib import Path
from task14_standard_kraken import read, write, parse_kraken


def details(run):
    start=time.time()
    dest=run/'profile-details'
    dest.mkdir()
    contigs={r['sequence_id']:r for r in read(run/'inputs/contig-map.tsv')}
    taxonomy={r['raw_taxid']:r for r in read(run/'summaries/taxonomy-calls.tsv')}
    approved={r['accession']:r for r in read(run/'inputs/accession-assignments.tsv')}
    profile_lookup={(r['accession'],r['input_mode'],r['confidence'],r['raw_taxid']):r for r in read(run/'summaries/taxid-profiles.tsv')}
    accounting={r['accession']:r for r in read(run/'inputs/sequence-accounting.tsv')}
    native=[]
    for conf in ['0','0.1']:
        with gzip.open(run/'raw'/f'native.confidence{conf}.output.gz','rt') as handle:
            for line in handle:
                acc,cid,taxid,length=parse_kraken(line,contigs,'native')
                p=profile_lookup[(acc,'native',conf,taxid)]
                native.append(dict(accession=acc,sequence_id=cid,confidence=conf,source_header=contigs[cid]['source_header'],length=length,ambiguous_bases=contigs[cid]['ambiguous_bases'],**taxonomy[taxid],framework_comparison=p['framework_comparison'],comparison_reason=p['comparison_reason'],project_species=approved[acc]['project_species'],project_subspecies=approved[acc]['project_subspecies'],approved_flags=approved[acc]['nomenclature_and_boundary_flags']))
    write(dest/'native-contig-calls.tsv',native,list(native[0]))
    # This map pass independently sums every represented source coordinate and
    # identifies all ambiguous fragments without holding 35.9m IDs in memory.
    ambiguous_ids={}
    map_count=map_bases=map_ambiguous_bases=0
    with gzip.open(run/'inputs/fragment-map.tsv.gz','rt') as handle:
        next(handle)
        for line in handle:
            fields=line.rstrip('\n').split('\t')
            map_count+=1
            map_bases+=int(fields[5])
            map_ambiguous_bases+=int(fields[7])
            if int(fields[7]):
                ambiguous_ids[fields[0]]=int(fields[7])
    assert map_count==sum(int(r['fragments']) for r in accounting.values())
    assert map_bases==sum(int(r['bases']) for r in accounting.values())
    assert map_ambiguous_bases==sum(int(r['ambiguous_bases']) for r in accounting.values())
    rows=[]
    for conf in ['0','0.1']:
        strata=collections.defaultdict(collections.Counter)
        with gzip.open(run/'raw'/f'fragments150.confidence{conf}.output.gz','rt') as handle:
            for line in handle:
                # The main summarizer already validates every ID, coordinate,
                # length and uniqueness; retain independent count arithmetic here.
                status,seqid,taxid,length,_=line.rstrip('\n').split('\t')
                length=int(length)
                acc=seqid.split('|',1)[0]
                length_stratum='full_150bp' if length==150 else 'tail_35_149bp' if length>=35 else 'tail_1_34bp'
                key=acc,length_stratum,int(seqid in ambiguous_ids)
                c=strata[key]
                c.update(fragments=1,bases=length,ambiguous_bases=ambiguous_ids.get(seqid,0),classified_fragments=int(status=='C'),unclassified_fragments=int(status=='U'),classified_bases=length*int(status=='C'),unclassified_bases=length*int(status=='U'))
                if taxid!='0':
                    c[taxonomy[taxid]['resolution']+'_fragments']+=1
        totals=collections.defaultdict(collections.Counter)
        for (acc,length_stratum,amb),c in sorted(strata.items()):
            totals[acc].update(fragments=c['fragments'],bases=c['bases'])
            if length_stratum=='tail_1_34bp':
                assert c['classified_fragments']==0, 'Unexpected call on a <k35 fragment'
            assert c['classified_fragments']+c['unclassified_fragments']==c['fragments']
            assert c['classified_bases']+c['unclassified_bases']==c['bases']
            rows.append(dict(accession=acc,confidence=conf,length_stratum=length_stratum,contains_ambiguous_bases=amb,**{key:c[key] for key in ['fragments','bases','ambiguous_bases','classified_fragments','unclassified_fragments','classified_bases','unclassified_bases','species_or_below_fragments','genus_only_fragments','below_genus_without_named_species_fragments','above_genus_or_unresolved_fragments']},stratum_unclassified_fragment_fraction=c['unclassified_fragments']/c['fragments'],stratum_unclassified_base_fraction=c['unclassified_bases']/c['bases'],assembly_fragment_denominator=accounting[acc]['fragments'],assembly_base_denominator=accounting[acc]['bases']))
        assert set(totals)==set(accounting)
        for acc,c in totals.items():
            assert c['fragments']==int(accounting[acc]['fragments']) and c['bases']==int(accounting[acc]['bases'])
    write(dest/'fragment-quality-strata.tsv',rows,list(rows[0]))
    pooled=collections.defaultdict(collections.Counter)
    assembly_fractions=collections.defaultdict(lambda: collections.defaultdict(float))
    for profile in profile_lookup.values():
        for dimension,category in [('resolution',profile['resolution']),('framework',profile['framework_comparison'])]:
            key=profile['input_mode'],profile['confidence'],dimension,category
            pooled[key].update(sequences=int(profile['sequences']),bases=int(profile['bases']))
            assembly_fractions[key][profile['accession']]+=float(profile['represented_base_fraction'])
    cohort=[]
    base_denominator=sum(int(a['bases']) for a in accounting.values())
    for (mode,conf,dimension,category),counter in sorted(pooled.items()):
        key=mode,conf,dimension,category
        values=[assembly_fractions[key].get(acc,0) for acc in accounting]
        sequence_denominator=sum(int(a['contigs' if mode=='native' else 'fragments']) for a in accounting.values())
        cohort.append(dict(input_mode=mode,confidence=conf,dimension=dimension,category=category,pooled_sequences=counter['sequences'],pooled_represented_bases=counter['bases'],pooled_sequence_denominator=sequence_denominator,pooled_base_denominator=base_denominator,pooled_sequence_fraction=counter['sequences']/sequence_denominator,pooled_represented_base_fraction=counter['bases']/base_denominator,assembly_denominator=len(accounting),assemblies_with_positive_profile_fraction=sum(v>0 for v in values),mean_assembly_base_fraction=statistics.mean(values),median_assembly_base_fraction=statistics.median(values),interpretation='sequence_call_profiles_not_genome_accuracy'))
    write(dest/'cohort-profile-summary.tsv',cohort,list(cohort[0]))
    (dest/'validation.json').write_text(json.dumps(dict(status='PASS',native_contig_call_rows=len(native),fragment_map_rows=map_count,fragment_map_bases=map_bases,ambiguous_bases=map_ambiguous_bases,ambiguous_fragment_ids=len(ambiguous_ids),all_shorter_than_k35_fragments_unclassified=True,elapsed_seconds=time.time()-start),indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',required=True,type=Path)
    a=p.parse_args()
    details(a.run)


if __name__=='__main__':
    main()
