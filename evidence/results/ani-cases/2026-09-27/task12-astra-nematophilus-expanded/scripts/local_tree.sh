#!/usr/bin/env bash
set -euo pipefail
cd ${PROJECT_ROOT}
RUN=results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded
REF="$1"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/python
IQ=${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3
/usr/bin/time -v -o "$RUN/logs/local-$REF.time" "$PY" "$RUN/scripts/local_alignment.py" --reference "$REF" --threads 4 > "$RUN/logs/local-$REF.stdout" 2> "$RUN/logs/local-$REF.stderr"
OUT="$RUN/local/$REF"
mkdir "$OUT/tree"
/usr/bin/time -v -o "$OUT/tree/primary.time" "$IQ" -s "$OUT/primary.core.fna" -m MFP --mset HKY,GTR --mrate G -T 4 -B 1000 --alrt 1000 --bnni --seed 20260927 --keep-ident --prefix "$OUT/tree/primary" > "$OUT/tree/primary.stdout" 2> "$OUT/tree/primary.stderr"
