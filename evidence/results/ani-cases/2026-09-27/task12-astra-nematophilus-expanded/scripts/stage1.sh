#!/usr/bin/env bash
set -euo pipefail
cd ${PROJECT_ROOT}
RUN=results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH="/sw/local/rocky9/noarch/rcc/software/gtdb-tk/gtdbtk_2.7.2/bin:$PATH"
export GTDBTK_DATA_PATH=/scratch/opendata/genomics/GTDB/releases/release232/
/usr/bin/time -v -o "$RUN/logs/identify.time" gtdbtk identify --batchfile "$RUN/inputs/marker.batch.tsv" --out_dir "$RUN/markers" --cpus 2 --write_single_copy_genes > "$RUN/logs/identify.stdout" 2> "$RUN/logs/identify.stderr"
SKANI=${USER_HOME}/sw/miniforge3/envs/LRE/bin/skani
"$SKANI" --version > "$RUN/logs/skani-version.txt"
"$SKANI" dist --help > "$RUN/logs/skani-help.txt"
/usr/bin/time -v -o "$RUN/logs/skani.time" "$SKANI" dist --ql "$RUN/inputs/new.list" --rl "$RUN/inputs/all.list" -t 2 -s 1 -o "$RUN/outputs/skani-forward.tsv" > "$RUN/logs/skani-forward.stdout" 2> "$RUN/logs/skani-forward.stderr"
/usr/bin/time -v -o "$RUN/logs/skani-reverse.time" "$SKANI" dist --ql "$RUN/inputs/all.list" --rl "$RUN/inputs/new.list" -t 2 -s 1 -o "$RUN/outputs/skani-reverse.tsv" > "$RUN/logs/skani-reverse.stdout" 2> "$RUN/logs/skani-reverse.stderr"
