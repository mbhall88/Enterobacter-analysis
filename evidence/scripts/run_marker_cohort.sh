#!/usr/bin/env bash
# GTDB-Tk identify + both align forms over any frozen cohort batchfile.
#
# Task 07's wrapper read one hardcoded path, `$RUN_DIR/inputs/pilot.abs.batch.tsv`. This one takes
# --batchfile, --run-dir and --threads, so the pilot and the full cohort run the same code with
# different arguments rather than the same code with the word "pilot" edited out of it.
#
# Both alignment forms are produced every time, and that is deliberate. SET-16 makes the untrimmed,
# occupancy-filtered alignment primary because GTDB's canonical mask cannot separate this genus --
# at pilot scale it left 183 distinct sequences among 388 genomes and made two pairs of differently
# named type anchors byte-identical. The masked alignment is still produced, because that ratio is
# a result to re-measure at full scale, not an assumption to carry over.
#
# Nothing here is on PATH by default; every path comes from the TESTED block of
# docs/analysis-execution-plan.md and is re-asserted.
#
# Usage: scripts/run_marker_cohort.sh --batchfile FILE --run-dir DIR [--threads N]
set -euo pipefail

PROJECT=${PROJECT_ROOT}
BATCH=""
RUN_DIR=""
THREADS="${SLURM_CPUS_PER_TASK:-4}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batchfile) BATCH="$2"; shift 2 ;;
    --run-dir)   RUN_DIR="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    *) echo "Error: unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -n "$BATCH"   ]] || { echo "Error: --batchfile is required" >&2; exit 2; }
[[ -n "$RUN_DIR" ]] || { echo "Error: --run-dir is required" >&2; exit 2; }
[[ -s "$BATCH"   ]] || { echo "Error: missing or empty batchfile $BATCH" >&2; exit 2; }

cd "$PROJECT"
module load gtdb-tk/2.7.2+232
export GTDBTK_DATA_PATH=/scratch/opendata/genomics/GTDB/releases/release232/
export PATH="/sw/local/rocky9/noarch/rcc/software/gtdb-tk/gtdbtk_2.7.2/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

LOGS="$RUN_DIR/logs"
mkdir -p "$LOGS"
for d in "$RUN_DIR/identify" "$RUN_DIR/canonical" "$RUN_DIR/full"; do
  [[ -e "$d" ]] && { echo "Error: refusing to overwrite existing $d" >&2; exit 2; }
done

{
  echo "host: $(hostname)"
  echo "slurm_job_id: ${SLURM_JOB_ID:-none}"
  echo "cpus: $THREADS"
  echo "gtdbtk: $(gtdbtk --version 2>&1 | head -1)"
  echo "GTDBTK_DATA_PATH: $GTDBTK_DATA_PATH"
  echo "batchfile: $BATCH"
  echo "batch_rows: $(wc -l < "$BATCH")"
  echo "batch_sha256: $(sha256sum "$BATCH" | cut -d' ' -f1)"
  echo "started: $(date -Is)"
} > "$LOGS/marker-environment.txt"

echo "== identify =="
/usr/bin/time -v -o "$LOGS/identify.time" gtdbtk identify \
  --batchfile "$BATCH" --out_dir "$RUN_DIR/identify" \
  --cpus "$THREADS" --write_single_copy_genes \
  > "$LOGS/identify.stdout" 2> "$LOGS/identify.stderr"

echo "== align: canonical (GTDB mask) — evidence only under SET-16 =="
/usr/bin/time -v -o "$LOGS/align-canonical.time" gtdbtk align \
  --identify_dir "$RUN_DIR/identify" --out_dir "$RUN_DIR/canonical" \
  --skip_gtdb_refs --cpus "$THREADS" \
  > "$LOGS/align-canonical.stdout" 2> "$LOGS/align-canonical.stderr"

echo "== align: untrimmed (--skip_trimming) — the primary form's input =="
/usr/bin/time -v -o "$LOGS/align-full.time" gtdbtk align \
  --identify_dir "$RUN_DIR/identify" --out_dir "$RUN_DIR/full" \
  --skip_gtdb_refs --skip_trimming --cpus "$THREADS" \
  > "$LOGS/align-full.stdout" 2> "$LOGS/align-full.stderr"

echo "finished: $(date -Is)" >> "$LOGS/marker-environment.txt"
echo "marker stage complete"
