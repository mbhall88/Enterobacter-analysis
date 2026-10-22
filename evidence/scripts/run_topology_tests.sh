#!/usr/bin/env bash
# Task 17 lines 3 and 5: do the data *reject* a competing resolution, or only fail to prefer it?
#
# IQ-TREE 3.1.1 is a TESTING INSTRUMENT here, not the tree method. SET-15 made FastTree the
# method for the marker tree and nothing produced by this script replaces that tree.
#
# For one alignment and a directory of IQ-TREE -g constraint trees (one per competing resolution,
# written by scripts/build_conflict_subsets.py):
#   1. the model: ModelFinder (--mset/--mrate restricted as stated below), or a given -m MODEL, or a
#      given partition scheme (-p) whose charpartition carries the models;
#   2. the unconstrained ML tree under that model;
#   3. the ML tree under each constraint, same model;
#   4. KH, SH, weighted KH/SH, ELW and AU on all of them together, 10,000 RELL replicates, with
#      model parameters estimated on the unconstrained ML tree (-t ML --tree-fix).
# In IQ-TREE 3.1.1 the brief's "-z ... -zb 10000 -au" is spelled --trees ... --test 10000
# --test-au; --test-weight adds the weighted KH/SH columns.
#
# --keep-ident everywhere: a constraint may name a sequence that IQ-TREE would otherwise collapse
# as identical, and every candidate must be scored on the identical taxon set.
#
# Usage: bash scripts/run_topology_tests.sh --alignment F --constraints DIR --out DIR
#          [--seqtype AA|DNA] [--model M | --partition F] [--threads N] [--seed S]
set -euo pipefail

PROJECT=${PROJECT_ROOT}
cd "$PROJECT"
IQTREE=${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

ALN="" CONS="" OUT="" MODEL="" PART="" SEQTYPE="AA"
THREADS="${SLURM_CPUS_PER_TASK:-4}" SEED=20260923
while [[ $# -gt 0 ]]; do
  case "$1" in
    --alignment) ALN="$2"; shift 2 ;;
    --constraints) CONS="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --partition) PART="$2"; shift 2 ;;
    --seqtype) SEQTYPE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -s "$ALN" && -d "$CONS" && -n "$OUT" ]] || { echo "need --alignment, --constraints, --out" >&2; exit 2; }
[[ -z "$MODEL" || -z "$PART" ]] || { echo "--model and --partition are exclusive" >&2; exit 2; }
if [[ -e "$OUT" ]]; then echo "refusing to overwrite $OUT" >&2; exit 2; fi
mkdir -p "$OUT/logs"

run() {  # run <label> <iqtree args...>
  local label="$1"; shift
  local start=$(date +%s)
  /usr/bin/time -v -o "$OUT/logs/$label.time" "$IQTREE" "$@" \
    > "$OUT/logs/$label.stdout" 2> "$OUT/logs/$label.stderr"
  echo -e "$label\t$(( $(date +%s) - start ))" >> "$OUT/stage-seconds.tsv"
}

common=(-s "$ALN" -st "$SEQTYPE" -T "$THREADS" --seed "$SEED" --keep-ident)
{
  echo "iqtree: $("$IQTREE" --version | head -1)"
  echo "host: $(hostname)"; echo "slurm_job_id: ${SLURM_JOB_ID:-none}"; echo "threads: $THREADS"
  echo "alignment: $ALN $(sha256sum "$ALN" | cut -d' ' -f1)"
  echo "seed: $SEED"; echo "started: $(date -Is)"
} > "$OUT/environment.txt"

if [[ -n "$PART" ]]; then
  modelargs=(-p "$PART")
  echo "model: partition scheme $PART $(sha256sum "$PART" | cut -d' ' -f1)" >> "$OUT/environment.txt"
  run ml "${common[@]}" "${modelargs[@]}" --prefix "$OUT/ml"
else
  if [[ -z "$MODEL" ]]; then
    # ModelFinder over the matrices that won on this project's marker data before (Task 07:
    # Q.INSECT, then Q.PLANT, LG) plus the GTDB/FastTree-family ones; free-rate and gamma both.
    if [[ "$SEQTYPE" == AA ]]; then
      mf=(-m MFP --mset LG,WAG,JTT,Q.PFAM,Q.INSECT,Q.PLANT,Q.YEAST --mrate G,I+G,R,I+R)
    else
      mf=(-m MFP --mset GTR,TIM2,TIM3,TVM,HKY,K3Pu --mrate G,I+G,R,I+R)
    fi
    run ml "${common[@]}" "${mf[@]}" --prefix "$OUT/ml"
    MODEL=$(sed -n 's/^Best-fit model according to BIC: //p' "$OUT/ml.iqtree" | head -1)
    [[ -n "$MODEL" ]] || { echo "no best-fit model in $OUT/ml.iqtree" >&2; exit 3; }
  else
    run ml "${common[@]}" -m "$MODEL" --prefix "$OUT/ml"
  fi
  modelargs=(-m "$MODEL")
  echo "model: $MODEL" >> "$OUT/environment.txt"
fi

cp "$OUT/ml.treefile" "$OUT/candidates.nwk"
printf 'ml\t%s\n' "$OUT/ml.treefile" > "$OUT/candidates.tsv"
for con in "$CONS"/*.nwk; do
  name=$(basename "$con" .nwk)
  run "con-$name" "${common[@]}" "${modelargs[@]}" -g "$con" --prefix "$OUT/con-$name"
  cat "$OUT/con-$name.treefile" >> "$OUT/candidates.nwk"
  printf '%s\t%s\n' "$name" "$OUT/con-$name.treefile" >> "$OUT/candidates.tsv"
done

run au "${common[@]}" "${modelargs[@]}" -t "$OUT/ml.treefile" --tree-fix \
  --trees "$OUT/candidates.nwk" --test 10000 --test-au --test-weight --prefix "$OUT/au"

echo "finished: $(date -Is)" >> "$OUT/environment.txt"
echo "topology tests complete: $OUT"
