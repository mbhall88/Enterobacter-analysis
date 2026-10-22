"""Run the reviewed postprocessing stages in frozen order after the public array completes."""
from common import RUN, ROOT, DESIGN, ART, EXPERT, PYTHON, read, write, sha, json_write, verify_object
from score_layers import logged_process


def stage(name, argv):
    log=RUN/'logs'/name
    log.mkdir(exist_ok=False)
    json_write(log/'command.json',dict(argv=argv,cwd=str(ROOT)))
    status=logged_process(argv,log)
    if status:raise RuntimeError(f'{name} exited {status}; no automatic retry')


def main():
    verify_object()
    for pin in read(RUN/'inputs/postprocessing-driver-pins.tsv')+read(RUN/'inputs/layer2-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('reviewed postprocessing input changed')
    stage('layer3-collection',[str(PYTHON),str(RUN/'source/collect.py'),'--layer','layer3'])
    stage('layer2-replay',[str(PYTHON),str(RUN/'source/layer2.py')])
    # Seal scorer inputs only once their output producers have completed.
    paths=[EXPERT,ART/'expert-sample-option-stratified-n100.tsv',
        DESIGN/'source/score_validation.py',DESIGN/'source/draw_disagreement_review.py',
        RUN/'inputs/postprocessing-driver-pins.tsv',RUN/'inputs/layer2-pins.tsv',
        RUN/'outputs/layer2-replay/summary.json']
    paths += [RUN/'outputs'/f'{layer}-collection'/name for layer in ['layer1','layer3']
              for name in ['attempt-ledger.tsv','summary.json','results.tsv']]
    paths += [RUN/'inputs'/f'{layer}-truth.tsv' for layer in ['layer1','layer3','layer3-held']]
    paths += [__file__,RUN/'source/score_layers.py',RUN/'source/common.py']
    write(RUN/'inputs/scoring-pins.tsv',[dict(path=str(p),sha256=sha(p)) for p in paths])
    stage('score-all-layers',[str(PYTHON),str(RUN/'source/score_layers.py')])
    json_write(RUN/'outputs/finalization-complete.json',dict(status='complete',source_sha256=sha(__file__),
        scoring_complete_sha256=sha(RUN/'outputs/scoring/complete.json'),production_assignments_enabled=False))

if __name__=='__main__':main()
