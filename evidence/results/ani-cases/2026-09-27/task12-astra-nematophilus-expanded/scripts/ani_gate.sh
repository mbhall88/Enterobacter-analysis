#!/usr/bin/env bash
set -euo pipefail
cd ${PROJECT_ROOT}
RUN=results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/tmp/nem-ani-audit
.work/envs/task12-pyani/bin/python scripts/audit_task12_anim_outputs.py --run "$RUN" --output "$RUN/anim/directional-audit" > "$RUN/logs/ani-audit.stdout" 2> "$RUN/logs/ani-audit.stderr"
python "$RUN/scripts/finish_ani.py" > "$RUN/logs/finish-ani.stdout" 2> "$RUN/logs/finish-ani.stderr"
