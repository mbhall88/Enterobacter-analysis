"""Adapt manifests to the frozen scorer; never reinterpret or relabel an outcome."""
import csv
import datetime
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
from common import RUN, ROOT, DESIGN, ART, EXPERT, PYTHON, read, write, json_write, sha, verify_object, load_module
csv.field_size_limit(1_000_000_000)
SCORER=DESIGN/'source/score_validation.py'


def logged_process(argv, log):
    timer=time.monotonic()
    record=dict(exit_status=None, started=datetime.datetime.now(datetime.timezone.utc).isoformat(), infrastructure_error='')
    try:
        with (log/'stdout.txt').open('x') as out,(log/'stderr.txt').open('x') as err:
            result=subprocess.run(argv,cwd=ROOT,stdout=out,stderr=err,check=False)
        record['exit_status']=result.returncode
        return result.returncode
    except Exception as error:
        record['infrastructure_error']=f'{type(error).__name__}: {error}'
        raise
    finally:
        record['seconds']=time.monotonic()-timer
        json_write(log/'finished.json',record)


def select_results(truth, ledger):
    by_id={r['genome_id']:r for r in ledger}
    if len(by_id)!=len(ledger):raise ValueError('duplicate result ledger identity')
    selected=[]
    for t in truth:
        r=by_id[t['genome_id']]
        if int(r['exit_status'])==0:
            if not r['report_tsv'] or not r['report_sha256']:raise ValueError('successful result missing report')
            selected.append({k:r[k] for k in ('genome_id','report_tsv','report_sha256')})
    return selected


def expert_name_breakdown(truth, ledger, metadata, scorer):
    """Supplementary grouping only; sealed labels are scored unchanged by the frozen function."""
    results={r['genome_id']:r for r in select_results(truth,ledger)}
    groups=defaultdict(list)
    for t in truth:
        r=results.get(t['genome_id'])
        report=scorer.load_result(r['report_tsv'],r['report_sha256']) if r else None
        groups[metadata[t['genome_id']]['ncbi_name_normalised']].append(scorer.score(t,report))
    return {name:scorer.summarise(rows) for name,rows in sorted(groups.items())}


def score_one(name, truth_path, ledger, groups):
    output=RUN/'outputs/scoring'; truth=read(truth_path)
    manifest=output/(name+'-results.tsv')
    write(manifest,select_results(truth,ledger),['genome_id','report_tsv','report_sha256'])
    argv=[str(PYTHON),str(SCORER),'--truth',str(truth_path),'--results',str(manifest),'--output',str(output/name)]
    for group in groups:argv+=['--group',group]
    log=RUN/'logs'/('score-'+name);log.mkdir(exist_ok=False)
    json_write(log/'command.json',dict(argv=argv,truth_sha256=sha(truth_path),scorer_sha256=sha(SCORER),
        results_manifest_sha256=sha(manifest),started=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    if logged_process(argv,log):raise RuntimeError('frozen scorer failed: '+name)
    # Author-facing per-genome table contains only coded IDs, never source names/accessions.
    codes={t['genome_id']:t['coded_id'] for t in truth}
    scored=read(output/name/'scored-genomes.tsv')
    safe=[]
    for r in scored:
        item={k:v for k,v in r.items() if k not in ['genome_id','block','component_id','sample_group','ncbi_name_normalised']}
        safe.append(dict(coded_id=codes[r['genome_id']],**item))
    write(output/name/'coded-outcomes.tsv',safe)
    return json.loads((output/name/'summary.json').read_text())


def main():
    verify_object()
    for pin in read(RUN/'inputs/scoring-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('scoring input/adapter drift')
    replay=json.loads((RUN/'outputs/layer2-replay/summary.json').read_text())
    if (replay['queries'],replay['settings'],replay['replays'],replay['saved_pairs'],replay['new_ANI']) != (1034,9,9306,143504,False):
        raise ValueError('Layer2 must complete before final scoring')
    if replay['production_assignments_enabled'] is not False:raise ValueError('unexpected positive assignment policy')
    for layer,count in [('layer1',946),('layer3',4857)]:
        audit=json.loads((RUN/'outputs'/f'{layer}-collection/summary.json').read_text())
        if not audit['all_accounted'] or audit['genomes']!=count:raise ValueError('incomplete native layer')
    output=RUN/'outputs/scoring';output.mkdir(exist_ok=False)
    ledger1=read(RUN/'outputs/layer1-collection/attempt-ledger.tsv')
    ledger3=read(RUN/'outputs/layer3-collection/attempt-ledger.tsv')
    near={r['genome_id']:r['near_clone'] for r in ledger1+ledger3}
    layer1=read(RUN/'inputs/layer1-truth.tsv')
    for r in layer1:r.update(near_clone=near[r['genome_id']],rule_design_exposed=str(bool(r['rule_design_role'])))
    path1=output/'layer1-truth.tsv';write(path1,layer1)
    summaries={}
    summaries['layer1']=score_one('layer1',path1,ledger1,[
        'truth_kind','partition','partition,truth_kind','label_status','saved_pairs_ani_class','truth_species_id','stratum',
        'assembly_level','reference_in_clone_block','reference_in_split_block','near_clone','rule_design_exposed',
        'partition,rule_design_exposed','partition,stratum','partition,label_status','partition,saved_pairs_ani_class'])
    for name in ['layer3','layer3-held']:
        truth=read(RUN/'inputs'/(name+'-truth.tsv'))
        for r in truth:r['near_clone']=near[r['genome_id']]
        path=output/(name+'-truth.tsv');write(path,truth)
        summaries[name]=score_one(name,path,ledger3,[
            'crosswalk_kind','qc_stratum','assembly_level','crosswalk_kind,qc_stratum,assembly_level',
            'expert_sample','near_clone','ncbi_name_normalised'])
    summaries['expert']=score_one('expert',EXPERT,ledger3,[
        'label_source','sample_group','qc_stratum','near_clone','subspecies_label_source'])
    scorer=load_module(SCORER,'frozen_expert_name_grouping')
    json_write(output/'expert-ncbi-name-breakdown.json',dict(expert_truth_sha256=sha(EXPERT),
        metadata_sha256=sha(RUN/'inputs/layer3-truth.tsv'),
        groups=expert_name_breakdown(read(EXPERT),ledger3,
            {t['genome_id']:t for t in read(RUN/'inputs/layer3-truth.tsv')},scorer)))
    # Sealed near-clone annotations stay immutable. A difference is documented, never rewritten.
    mismatches=[dict(coded_id=t['coded_id'],sealed_near_clone=t['near_clone'],native_near_clone=near[t['genome_id']])
                for t in read(EXPERT) if t['near_clone']!=near[t['genome_id']]]
    write(output/'expert-near-clone-check.tsv',mismatches,['coded_id','sealed_near_clone','native_near_clone'])
    draw=DESIGN/'source/draw_disagreement_review.py'
    argv=[str(PYTHON),str(draw),'--scored',str(output/'layer3/scored-genomes.tsv'),
        '--truth',str(output/'layer3-truth.tsv'),'--expert-sample',str(ART/'expert-sample-option-stratified-n100.tsv'),
        '--output',str(output/'disagreement-draw')]
    log=RUN/'logs/draw-disagreements';log.mkdir(exist_ok=False)
    json_write(log/'command.json',dict(argv=argv,script_sha256=sha(draw)))
    if logged_process(argv,log):raise RuntimeError('frozen disagreement draw failed')
    json_write(output/'complete.json',dict(layers=list(summaries),scorer_sha256=sha(SCORER),expert_truth_sha256=sha(EXPERT),
        expert_near_clone_mismatches=len(mismatches),positive_assignments_enabled=False,script_sha256=sha(__file__)))
    print(json.dumps({name:summary['groups'][next(k for k in summary['groups'] if '|' not in k)]['species_outcomes']
                      for name,summary in summaries.items()}))

if __name__=='__main__':main()
