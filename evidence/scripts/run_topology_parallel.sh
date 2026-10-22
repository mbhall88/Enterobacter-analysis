#!/usr/bin/env bash
# Task 17 line 5: the same topology test as run_topology_tests.sh, split so each tree search is
# its own Slurm job. A 3-Mb core concatenation makes one search take hours, and the unconstrained
# search and one search per competing resolution are independent.
#
# IQ-TREE 3.1.1 is a TESTING INSTRUMENT here; SET-15's FastTree tree is not replaced.
#
#   --mode search  one ML search: unconstrained (no --constraint) -> DIR/ml.*, or constrained to
#                  one IQ-TREE -g tree -> DIR/con-<name>.*
#   --mode au      every DIR/ml.treefile and DIR/con-*.treefile together: KH, SH, weighted KH/SH,
#                  ELW, AU, 10,000 RELL replicates, model parameters estimated on the unconstrained
#                  ML tree (-t ml --tree-fix); writes DIR/candidates.tsv in the form
#                  scripts/summarise_topology_tests.py reads.
# The model comes only from --partition (a charpartition with a model per charset), so every search
# and the test use the identical model.
#
# Usage: bash scripts/run_topology_parallel.sh --mode search|au --alignment F --partition F
#          --dir DIR [--constraint F] [--seqtype DNA|AA] [--threads N] [--seed S]
set -euo pipefail
PROJECT=${PROJECT_ROOT}
cd "$PROJECT"
IQTREE=${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

MODE="" ALN="" PART="" DIR="" CON="" SEQTYPE="DNA" SEED=20260923
THREADS="${SLURM_CPUS_PER_TASK:-4}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --alignment) ALN="$2"; shift 2 ;;
    --partition) PART="$2"; shift 2 ;;
    --dir) DIR="$2"; shift 2 ;;
    --constraint) CON="$2"; shift 2 ;;
    --seqtype) SEQTYPE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -s "$ALN" && -s "$PART" && -n "$DIR" && "$MODE" =~ ^(search|au)$ ]] \
  || { echo "need --mode search|au, --alignment, --partition, --dir" >&2; exit 2; }
mkdir -p "$DIR/logs"
common=(-s "$ALN" -st "$SEQTYPE" -p "$PART" -T "$THREADS" --seed "$SEED" --keep-ident)

run() {
  local label="$1"; shift
  local start=$(date +%s)
  /usr/bin/time -v -o "$DIR/logs/$label.time" "$IQTREE" "$@" \
    > "$DIR/logs/$label.stdout" 2> "$DIR/logs/$label.stderr"
  echo -e "$label\t$(( $(date +%s) - start ))\t${SLURM_JOB_ID:-none}" >> "$DIR/stage-seconds.tsv"
}

if [[ "$MODE" == search ]]; then
  if [[ -z "$CON" ]]; then
    prefix="$DIR/ml" label=ml
  else
    name=$(basename "$CON" .nwk); prefix="$DIR/con-$name" label="con-$name"
  fi
  [[ -e "$prefix.treefile" ]] && { echo "refusing to overwrite $prefix.treefile" >&2; exit 2; }
  {
    echo "iqtree: $("$IQTREE" --version | head -1)"; echo "alignment: $ALN $(sha256sum "$ALN" | cut -d' ' -f1)"
    echo "partition: $PART $(sha256sum "$PART" | cut -d' ' -f1)"; echo "seed: $SEED"
  } > "$DIR/environment.txt.tmp.$label"
  if [[ -z "$CON" ]]; then
    run "$label" "${common[@]}" --prefix "$prefix"
  else
    run "$label" "${common[@]}" -g "$CON" --prefix "$prefix"
  fi
  mv "$DIR/environment.txt.tmp.$label" "$DIR/environment-$label.txt"
  echo "search complete: $prefix.treefile"
  exit 0
fi

[[ -s "$DIR/ml.treefile" ]] || { echo "no unconstrained tree in $DIR" >&2; exit 3; }
[[ -e "$DIR/au.iqtree" ]] && { echo "refusing to overwrite $DIR/au.iqtree" >&2; exit 2; }
cp "$DIR/ml.treefile" "$DIR/candidates.nwk"
printf 'ml\t%s\n' "$DIR/ml.treefile" > "$DIR/candidates.tsv"
for f in "$DIR"/con-*.treefile; do
  name=$(basename "$f" .treefile); name="${name#con-}"
  cat "$f" >> "$DIR/candidates.nwk"
  printf '%s\t%s\n' "$name" "$f" >> "$DIR/candidates.tsv"
done
run au "${common[@]}" -t "$DIR/ml.treefile" --tree-fix --trees "$DIR/candidates.nwk" \
  --test 10000 --test-au --test-weight --prefix "$DIR/au"
{
  echo "iqtree: $("$IQTREE" --version | head -1)"; echo "alignment: $ALN $(sha256sum "$ALN" | cut -d' ' -f1)"
  echo "model: partition $PART"; echo "candidates: $(wc -l < "$DIR/candidates.tsv")"
  echo "finished: $(date -Is)"
} > "$DIR/environment.txt"
echo "topology test complete: $DIR"
