"""Reuse the byte-pinned original worker for only documented incomplete inputs."""
from pathlib import Path
import sys
import json
ORIGINAL=Path('${PROJECT_ROOT}/results/tool-planning/2026-10-07/validation-runs-astra')
RECOVERY=ORIGINAL.with_name('validation-quota-recovery-astra')
sys.path.insert(0,str(ORIGINAL/'source'))
import run_batch
from common import read,sha

if __name__=='__main__':
    for pin in read(RECOVERY/'inputs/recovery-driver-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('recovery wrapper changed')
    for pin in read(RECOVERY/'inputs/execution-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('recovery input changed')
    original={r['genome_id']:r for r in read(ORIGINAL/'inputs/execution-manifest.tsv')}
    successes={r['genome_id'] for r in read(ORIGINAL/'outputs/prearchive-successful-ledger.tsv')}
    todo=read(RECOVERY/'inputs/execution-manifest.tsv')
    reason_rows=read(RECOVERY/'inputs/recovery-dispositions.tsv')
    reasons={r['genome_id']:r for r in reason_rows}
    if (len({r['genome_id'] for r in todo})!=len(todo) or len(reasons)!=len(reason_rows)
        or set(reasons)!={r['genome_id'] for r in todo}):
        raise ValueError('recovery dispositions must match unique inputs exactly')
    if {r['genome_id'] for r in todo}!={k for k,r in original.items() if r['layer']=='layer3'}-successes:
        raise ValueError('recovery set must exactly cover unsuccessful/unattempted public inputs')
    for row in todo:
        ident=row['genome_id']
        if row!=original[ident] or reasons[ident]['previous_state'] not in {'unattempted','interrupted_by_quota_shutdown','failed_disk_quota'}:
            raise ValueError('recovery input/reason differs from authorised infrastructure scope')
        previous=ORIGINAL/'queries'/row['coded_id']/'attempt01/finished.json'
        if previous.exists() and json.loads(previous.read_text()).get('exit_status')==0:
            raise ValueError('refusing to rerun successful original input')
        native=previous.parent/'analyze/run.json'
        if native.exists() and json.loads(native.read_text()).get('status')=='complete':
            raise ValueError('completed native output requires recovery, not rerunning')
    run_batch.RUN=RECOVERY
    run_batch.main()
