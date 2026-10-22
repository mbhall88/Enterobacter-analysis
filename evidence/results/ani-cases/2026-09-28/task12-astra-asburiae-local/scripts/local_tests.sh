#!/usr/bin/env bash
set -euo pipefail
cd ${PROJECT_ROOT}
RUN=results/ani-cases/2026-09-28/task12-astra-asburiae-local
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/python "$RUN/scripts/local_tests.py" --reference "$1" > "$RUN/logs/tests-$1.stdout" 2> "$RUN/logs/tests-$1.stderr"
