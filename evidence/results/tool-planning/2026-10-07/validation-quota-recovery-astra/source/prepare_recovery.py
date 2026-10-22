"""Record the quota incident and freeze the exact incomplete-input resume set."""
import sys,json,datetime,shutil
from pathlib import Path
ROOT=Path('${PROJECT_ROOT}')
ORIGINAL=ROOT/'results/tool-planning/2026-10-07/validation-runs-astra'
RECOVERY=ORIGINAL.with_name('validation-quota-recovery-astra')
LOCAL=Path('/tmp/enterobacter-validation-quota-recovery-2026-10-07')
sys.path.insert(0,str(ORIGINAL/'source'))
from common import read,write,sha,json_write,verify_object


def copy_new(source,destination):
    with source.open('rb') as inp,destination.open('xb') as out:shutil.copyfileobj(inp,out)
    if sha(source)!=sha(destination):raise ValueError('copy differs')


def main():
    verify_object();archive=json.loads((ORIGINAL/'outputs/ancillary-archive.json').read_text())
    if archive['status']!='complete':raise ValueError('archival recovery incomplete')
    good={r['genome_id'] for r in read(ORIGINAL/'outputs/prearchive-successful-ledger.tsv')}
    original=read(ORIGINAL/'inputs/execution-manifest.tsv');todo=[];reasons=[]
    for r in original:
        if r['genome_id'] in good:continue
        if r['layer']!='layer3':raise ValueError('previously completed Layer1 missing')
        directory=ORIGINAL/'queries'/r['coded_id']/'attempt01';finished=directory/'finished.json'
        native=directory/'analyze/run.json';started=directory/'started.json';stderr=directory/'stderr.txt'
        if native.exists() and json.loads(native.read_text()).get('status')=='complete':raise ValueError('recover native success without rerunning')
        if finished.exists():
            result=json.loads(finished.read_text())
            if result['exit_status']==0:raise ValueError('successful result not in archived ledger')
            if not stderr.exists() or 'Disk quota exceeded' not in stderr.read_text():raise ValueError('non-quota failure needs specific adjudication')
            state='failed_disk_quota';why='native exit2 with explicit Errno122; identical-input infrastructure retry'
        elif started.exists():state='interrupted_by_quota_shutdown';why='public array cancelled on quota incident; native completion absent; identical-input retry'
        elif not directory.exists():state='unattempted';why='original input never started'
        else:raise ValueError('unexpected reserved output without attempt record')
        todo.append(r);reasons.append(dict(genome_id=r['genome_id'],coded_id=r['coded_id'],previous_state=state,reason=why,
            original_attempt_directory=str(directory),started_sha256=sha(started) if started.exists() else '',
            finished_sha256=sha(finished) if finished.exists() else '',stderr_sha256=sha(stderr) if stderr.exists() else ''))
    if len(good)!=5045 or len(todo)!=758:raise ValueError('unexpected recovery counts')
    from collections import Counter
    counts=dict(Counter(r['previous_state'] for r in reasons))
    if counts!={'unattempted':740,'interrupted_by_quota_shutdown':16,'failed_disk_quota':2}:raise ValueError('unexpected recovery causes')
    RECOVERY.mkdir(exist_ok=False)
    for name in ['source','inputs','logs','outputs','reviews']:(RECOVERY/name).mkdir()
    claim=ROOT/'docs/task-status/15-validation-quota-recovery-astra.md'
    with claim.open('x') as handle:handle.write('# Task15 infrastructure recovery — file quota\n\nOwner: Codex Astra (`gpt-6-astra`), 7 October2026.\nStatus: running; frozen validation step4 continuation under design section9.\n\nRun: `results/tool-planning/2026-10-07/validation-quota-recovery-astra/`.\nOriginal: `results/tool-planning/2026-10-07/validation-runs-astra/`.\n\nThe one-million-file hard quota caused2explicit failures and16interrupted attempts;740public inputs had not started. The5045successful inputs are reused. Completed ancillary bytes are retained in a fully verified archive; key reports/evidence remain loose. No classifier, policy, reference, query or sealed label changed. Resume only the758documented missing inputs with the same frozen command. No RDM access, commits, other-agent changes or successful-input reruns. Completion requires all original5803inputs accounted for, archive-aware integrity verification, frozen Layer2/scoring and updated Task15 plan row.\n')
    for name in ['archive_completed.py','test_archive_fixture.py','recovery_worker.py','collect_recovery.py','finalize_recovery.py','prepare_recovery.py']:
        copy_new(LOCAL/name,RECOVERY/'source'/name)
    for name in ['archive.stdout','archive.stderr','archive.resources','archive-result.json','archive-verified-checkpoint.json',
                 'synthetic-archive-test.json','synthetic-archive-test.stdout','synthetic-archive-test.stderr','attempt-inventory.json','interrupted-native-status.json']:
        copy_new(LOCAL/name,RECOVERY/'logs'/name)
    copy_new(Path('/tmp/enterobacter-validation-quota-incident-2026-10-07.json'),RECOVERY/'logs/quota-incident.json')
    write(RECOVERY/'inputs/execution-manifest.tsv',todo)
    write(RECOVERY/'inputs/recovery-dispositions.tsv',reasons)
    paths=[ORIGINAL/'inputs/execution-manifest.tsv',ORIGINAL/'inputs/execution-pins.tsv',ORIGINAL/'inputs/driver-pins.tsv',
        ORIGINAL/'outputs/prearchive-successful-ledger.tsv',ORIGINAL/'outputs/ancillary-archive.json',
        ORIGINAL/'outputs/archive-verified-checkpoint.json',ORIGINAL/'outputs/ancillary-inventory.tsv',
        RECOVERY/'inputs/execution-manifest.tsv',RECOVERY/'inputs/recovery-dispositions.tsv']
    paths += [Path(r['path']) for r in read(ORIGINAL/'inputs/execution-pins.tsv')]
    write(RECOVERY/'inputs/execution-pins.tsv',[dict(path=str(p),sha256=sha(p)) for p in dict.fromkeys(paths)])
    paths=[ORIGINAL/'source'/n for n in ['common.py','run_batch.py','collect.py','finalize.py','score_layers.py','layer2.py']]
    paths += sorted((RECOVERY/'source').glob('*.py'))
    paths += [RECOVERY/'inputs/execution-pins.tsv']
    pins=[dict(path=str(p),sha256=sha(p)) for p in paths]
    write(RECOVERY/'inputs/recovery-driver-pins.tsv',pins);write(RECOVERY/'inputs/driver-pins.tsv',pins)
    json_write(RECOVERY/'inputs/preparation.json',dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),inputs=758,
        previous_states=counts,successful_inputs_reused=5045,successful_inputs_rerun=0,
        original_manifest_sha256=sha(ORIGINAL/'inputs/execution-manifest.tsv'),recovery_manifest_sha256=sha(RECOVERY/'inputs/execution-manifest.tsv'),
        archive_sha256=archive['archive_sha256'],frozen_classifier_changed=False,rdm_touched=False))
    (RECOVERY/'README.md').write_text('# Task15 quota recovery\n\nCodex Astra (`gpt-6-astra`), 7October2026.\n\nStatus: prepared to resume758public inputs:740unattempted,16cancel-interrupted,2explicit disk-quota failures. All5045successful results are reused, and all earlier attempts retained. The frozen validation command and scientific inputs are unchanged.\n\nCompleted ancillary files are retained at `../validation-runs-astra/outputs/completed-native-ancillary.tar.gz`; every member is checked against native hashes. Seven core native outputs and all outer attempt logs remain at their original paths. `ancillary-inventory.tsv` records exact archive paths and hashes. To inspect ancillary files, extract selected members into temporary space; restoring the entire archive on this fileset could exceed its file quota again.\n\nThe recovery worker uses the original reviewed code with a new output root. It verifies exact original query rows and documented infrastructure reasons and refuses successful original/native-complete inputs before running. The merged collector verifies archive plus loose files, enforces the unchanged4857public denominator and retains attempt provenance. The original reviewed finalizer then runs Layer2 and scoring. Scheduler IDs and completion results will follow.\n')
    print(json.dumps(dict(run=str(RECOVERY),inputs=len(todo),previous_states=counts)))

if __name__=='__main__':main()
