#!/usr/bin/env bash
# Task 11 -- the three-method comparison at full scale, under SET-18.
#
# Everything panel-specific lives in config/tree-comparison-full.json, which is derived from the
# Task 10 configuration by commands/make_full_comparison_config.py. That is the point of the
# design: "full" is a configuration change, not a text substitution, and the comparison code is
# byte-identical to the code whose split machinery reproduced Task 08's and Task 09's own figures.
set -euo pipefail
PROJECT=${PROJECT_ROOT}
RUN="$PROJECT/results/full-cohort/2026-09-22/task11-opus-full"
cd "$PROJECT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

/usr/bin/time -v -o "$RUN/logs/comparison.time" python3 scripts/compare_phylogenies.py \
  --config config/tree-comparison-full.json \
  --output "$RUN/comparison" --root . \
  > "$RUN/logs/comparison.stdout" 2> "$RUN/logs/comparison.stderr"

echo "comparison written to $RUN/comparison"
