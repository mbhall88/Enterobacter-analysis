#!/usr/bin/env bash
set -euo pipefail
cd ${PROJECT_ROOT}
TASK14_RUN=${PROJECT_ROOT}/results/standard-kraken/2026-10-03/task14-sol61-whole-fragments
TASK14_PYTHON=${USER_HOME}/sw/miniforge3/bin/python
trap 'TASK14_EXIT=$?; printf "%s\n" "$TASK14_EXIT" > "$TASK14_RUN/logs/worker.exit"' EXIT
/usr/bin/time -v -o "$TASK14_RUN/logs/worker.time" "$TASK14_PYTHON" "$TASK14_RUN/source-snapshot/task14_standard_kraken.py" pipeline --output "$TASK14_RUN" --threads 16
