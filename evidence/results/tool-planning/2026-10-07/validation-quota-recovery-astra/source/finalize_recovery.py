"""Use the reviewed finalizer with an archive-aware merged public collector."""
from pathlib import Path
import sys
ORIGINAL=Path('${PROJECT_ROOT}/results/tool-planning/2026-10-07/validation-runs-astra')
RECOVERY=ORIGINAL.with_name('validation-quota-recovery-astra')
sys.path.insert(0,str(ORIGINAL/'source'))
import finalize
from common import read,sha,PYTHON

if __name__=='__main__':
    for pin in read(RECOVERY/'inputs/recovery-driver-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('recovery wrapper changed')
    original_stage=finalize.stage
    def stage(name,argv):
        if name=='layer3-collection':argv=[str(PYTHON),str(RECOVERY/'source/collect_recovery.py')]
        return original_stage(name,argv)
    finalize.stage=stage
    finalize.main()
