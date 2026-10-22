"""Prespecified development replay of saved ANI, without new native analyses or tuning."""
import copy
import gzip
import hashlib
import json
import time
from collections import Counter, defaultdict
from itertools import product
from common import ROOT, RUN, DESIGN, PINS, read, write, sha, json_write, load_module, verify_object
from prepare import coded, truth_base

CAL = ROOT / 'results/tool-planning/2026-10-06/rank-calibration-astra'
BASE = ROOT / 'results/tool-planning/2026-10-06/approved-scope-followthrough-astra/development-baseline'
PAIR_SHA = '8f003c266a85d09ad8814718391890cbed41674b4a4bd2a723568f0d6cc590ba'
QUERY_SHA = 'c93888e29f1b976f933bb1d2af1a8ccaa57560113efd012d1f5d38af4de57544'


def masked_context(context, accessions):
    """Run-local subset; the package, context and immutable source files stay untouched."""
    subset = copy.copy(context)
    subset._references = {k: r for k, r in context._references.items() if r['accession'] in accessions}
    return subset


def replay_row(hits, query, refs, context, rank, config):
    from enterobacter_classifier.evidence import collect_evidence
    from enterobacter_classifier.eligibility import evaluate_eligibility, primary_inventory
    from enterobacter_classifier.reporting import _eligibility_columns
    inventory = primary_inventory(context)
    refs_by_accession = {r['accession']: r for r in refs}
    rows = context.evidence_references()
    # Frozen package IDs are retained; historical paths explicitly map saved pairs to those IDs.
    for r in rows:
        ref = context._references[r['reference_id']]
        r['search_path'] = refs_by_accession[ref['accession']]['whole_fasta']
        if ref['source']['sha256'] != refs_by_accession[ref['accession']]['whole_sha256']:
            raise ValueError('saved reference bytes differ from package')
    provenance = dict(input_kind='assembly', sequence_scope='whole_assembly_as_supplied', engine='skani',
        engine_version='0.3.2', search_exit_code=0, query_sha256=query['whole_sha256'],
        search_parameters_record=str(BASE / 'design.json'), raw_hits_sha256=PAIR_SHA,
        reference_package_sha256=PINS['package'][1], reference_component_id='primary_whole',
        execution_provenance='saved_development_pairs_prespecified_inventory_exclusion_replay')
    evidence = collect_evidence(hits, rows, query_search_path=query['whole_fasta'], provenance=provenance)
    digest = hashlib.sha256(json.dumps(evidence, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    result = evaluate_eligibility(evidence, rank, config, reference_inventory=inventory, evidence_sha256=digest)
    report = _eligibility_columns(dict(experimental_eligibility=result, assignment=evidence['assignment'],
        nearest_candidates=[h for h in evidence['ranked_reference_evidence'] if h['reference_id'] in evidence['nearest_reference_ids']]))
    report['query_sha256'] = query['whole_sha256']
    return report, result, digest


def main():
    started=time.monotonic()
    verify_object()
    for layer,count in [('layer1',946),('layer3',4857)]:
        audit=json.loads((RUN/'outputs'/f'{layer}-collection/summary.json').read_text())
        if not audit['all_accounted'] or audit['genomes'] != count:
            raise ValueError('native validation layers must finish before Layer2')
    for r in read(RUN/'inputs/layer2-pins.tsv'):
        if sha(r['path']) != r['sha256']:raise ValueError('Layer2 replay input changed')
    if sha(BASE/'hits.tsv') != PAIR_SHA or sha(CAL/'inputs/queries.tsv') != QUERY_SHA:
        raise ValueError('saved development input digest mismatch')
    old=load_module(CAL/'source/analyze.py','frozen_calibration_exclusions')
    scorer=load_module(DESIGN/'source/score_validation.py','frozen_layer2_scorer')
    from enterobacter_classifier.reporting import ReportingContext
    context=ReportingContext(PINS['package'][0],expected_sha256=PINS['package'][1])
    species={c['display_name']:c['concept_id'] for c in context._concepts.values() if c['rank']=='species'}
    subspecies={c['display_name']:c['concept_id'] for c in context._concepts.values() if c['rank']=='subspecies'}
    queries,refs=read(CAL/'inputs/queries.tsv'),read(CAL/'inputs/references.tsv')
    if len(queries)!=1034 or len(refs)!=140:raise ValueError('wrong development dimensions')
    if {r['accession'] for r in refs}!={r['accession'] for r in context._references.values()}:
        raise ValueError('saved and package reference identities differ')
    if any(r['partition']!='development' for r in queries+refs):raise ValueError('reserve in development replay')
    raw=defaultdict(list);qp={q['whole_fasta']:q for q in queries};rp={r['whole_fasta']:r for r in refs};pairs=set()
    for h in read(BASE/'hits.tsv'):
        q,r=qp[h['Query_file']],rp[h['Ref_file']]
        pair=q['accession'],r['accession']
        if pair in pairs:raise ValueError('duplicate saved pair')
        pairs.add(pair);raw[q['accession']].append(h)
    if len(pairs)!=143504:raise ValueError('saved pair count differs')
    design=json.loads((CAL/'design.json').read_text())
    if (design['modes']!=['self_alias','clone_block','split_block']
        or design['depletion_arms']!=['represented','omit_reviewed_species','omit_reviewed_subspecies']):
        raise ValueError('prespecified settings changed')
    rank=json.loads(PINS['rank'][0].read_text());config=json.loads(PINS['config'][0].read_text())
    output=RUN/'outputs/layer2-replay';output.mkdir(exist_ok=False)
    summaries=defaultdict(list);rows=[];excluded=[]
    with gzip.open(output/'eligibility-records.jsonl.gz','xt') as handle:
        for mode,arm in product(design['modes'],design['depletion_arms']):
            for q in queries:
                removed={r['accession']:old.exclude(q,r,mode,arm) for r in refs if old.exclude(q,r,mode,arm)}
                remaining=[r for r in refs if r['accession'] not in removed]
                subset=masked_context(context,{r['accession'] for r in remaining})
                hits=[h for h in raw[q['accession']] if rp[h['Ref_file']]['accession'] not in removed]
                report,result,digest=replay_row(hits,q,remaining,subset,rank,config)
                truth=truth_base(q['accession'],'layer2',q['whole_sha256'],q['clone_block'])
                named=q['named_target']=='1'
                truth.update(truth_kind='named' if named else 'context',label_source='project_label',
                    truth_species_id=species[q['project_species']] if named else '',
                    truth_subspecies_id=subspecies[q['project_species']+' subsp. '+q['project_subspecies']] if named and q['project_subspecies'] else '')
                outcome=scorer.score(truth,report)
                applicable=arm=='represented' or (named and (arm=='omit_reviewed_species' or bool(q['project_subspecies'])))
                info=dict(mode=mode,arm=arm,omission_applicable=str(applicable),coded_id=coded(q['accession'],'l2'),
                    project_species=q['project_species'],project_subspecies=q['project_subspecies'],qc_stratum=old.quality(q),
                    references_available=len(remaining),reported_hits=len(hits),excluded_reference_count=len(removed),
                    subspecies_reason=result['rank_evidence']['subspecies']['reason'])
                rows.append(dict(info,**outcome))
                # Inapplicable omission cells are retained as such, never pooled into arm summaries.
                if applicable:
                    for label in ['all','quality='+info['qc_stratum'],'species='+info['project_species'],
                                  'species='+info['project_species']+'|quality='+info['qc_stratum']]:
                        summaries[mode+'|'+arm+'|'+label].append(outcome)
                excluded.extend(dict(mode=mode,arm=arm,coded_id=info['coded_id'],reference=r,reason=why) for r,why in removed.items())
                handle.write(json.dumps(dict(mode=mode,arm=arm,coded_id=info['coded_id'],evidence_digest=digest,
                    report=report,eligibility=result),sort_keys=True)+'\n')
    write(output/'scored-development.tsv',rows);write(output/'reference-exclusions.tsv',excluded)
    write(output/'coded-development-outcomes.tsv', [{k:v for k,v in r.items() if k!='genome_id'} for r in rows])
    if len(rows)!=9306:raise ValueError('development replay incomplete')
    json_write(output/'summary.json',dict(queries=1034,settings=9,replays=len(rows),saved_pairs=143504,
        new_ANI=False,production_assignments_enabled=False,interpretation='development outcomes, not independent accuracy',
        scorer_sha256=sha(DESIGN/'source/score_validation.py'),seconds=time.monotonic()-started,
        groups={k:scorer.summarise(v) for k,v in sorted(summaries.items())}))
    print(json.dumps(dict(replays=len(rows),saved_pairs=143504,new_ANI=False,seconds=time.monotonic()-started)))

if __name__=='__main__':main()
