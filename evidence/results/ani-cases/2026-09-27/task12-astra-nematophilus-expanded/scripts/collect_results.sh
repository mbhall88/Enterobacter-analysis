#!/usr/bin/env bash
set -euo pipefail
cd ${PROJECT_ROOT}
RUN=results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$RUN/logs/matplotlib-cache"
PY=${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/python
"$PY" "$RUN/scripts/validate_masked_fixture.py" > "$RUN/logs/final-fixture.stdout" 2> "$RUN/logs/final-fixture.stderr"
"$PY" "$RUN/scripts/summarize_trees.py" > "$RUN/logs/tree-summary.stdout" 2> "$RUN/logs/tree-summary.stderr"
for condition in primary masked RY; do
    "$PY" "$RUN/scripts/plot_local_trees.py" --condition "$condition" > "$RUN/logs/plot-$condition.stdout" 2> "$RUN/logs/plot-$condition.stderr"
done
"$PY" "$RUN/scripts/validate_final.py" > "$RUN/logs/final-validation.stdout" 2> "$RUN/logs/final-validation.stderr"
