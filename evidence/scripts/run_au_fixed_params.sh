#!/usr/bin/env bash
# Task 17: AU/KH/SH test with the substitution model held at one fitted optimum.
#
# run_topology_parallel.sh --mode au re-estimates every partition's GTR+F+R4 parameters from
# scratch on the unconstrained tree before scoring the candidates.  On the C3 panel that fresh
# estimate settled 3,670 lnL below the search's own optimum for the same tree (tree length 0.0821
# against 0.0854), and the candidates' ranking under it differed from their search scores.  This
# wrapper scores every candidate under ONE parameter set: that of the candidate whose own search
# reached the highest log-likelihood (read from its .iqtree, chosen before any test is run).
# IQ-TREE fixes parameters written inside braces in the model string, so passing that search's
# .best_model.nex as the partition file leaves branch lengths (and partition rates) free and
# nothing else.
#
#   run_au_fixed_params.sh --alignment ALN --dir SEARCH_DIR --out OUT \
#       --candidate NAME=TREEFILE [--candidate ...] [--threads N]
#
# Each TREEFILE must be a finished search in SEARCH_DIR (its .iqtree and .best_model.nex beside
# it).  Writes OUT/candidates.tsv, OUT/candidates.nwk, OUT/params-from.txt, OUT/au.* and
# OUT/environment.txt (with a `finished:` line) so summarise_topology_tests.py can read OUT.
set -euo pipefail

IQTREE=${IQTREE:-${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3}
ALN="" DIR="" OUT="" THREADS=16 SEED=20260923
CANDS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --alignment) ALN=$2; shift 2 ;;
    --dir) DIR=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --candidate) CANDS+=("$2"); shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$ALN" && -n "$DIR" && -n "$OUT" && ${#CANDS[@]} -ge 2 ]] || { echo "missing arguments" >&2; exit 2; }
[[ -e "$OUT/au.iqtree" ]] && { echo "refusing to overwrite $OUT" >&2; exit 1; }
mkdir -p "$OUT"

: > "$OUT/candidates.tsv"; : > "$OUT/candidates.nwk"; : > "$OUT/search-lnl.tsv"
best_name="" best_lnl="" best_prefix=""
for c in "${CANDS[@]}"; do
  name=${c%%=*} tree=${c#*=}
  prefix=${tree%.treefile}
  [[ -s "$tree" && -s "$prefix.iqtree" && -s "$prefix.best_model.nex" ]] \
    || { echo "candidate $name: $tree is not a finished search" >&2; exit 1; }
  lnl=$(sed -n 's/^Log-likelihood of the tree: \(-[0-9.]*\).*/\1/p' "$prefix.iqtree")
  [[ -n "$lnl" ]] || { echo "no log-likelihood in $prefix.iqtree" >&2; exit 1; }
  printf '%s\t%s\n' "$name" "$tree" >> "$OUT/candidates.tsv"
  cat "$tree" >> "$OUT/candidates.nwk"
  if [[ -z "$best_lnl" ]] || awk -v a="$lnl" -v b="$best_lnl" 'BEGIN{exit !(a>b)}'; then
    best_name=$name best_lnl=$lnl best_prefix=$prefix
  fi
  printf '%s\t%s\n' "$name" "$lnl" >> "$OUT/search-lnl.tsv"
done
printf 'parameters_from\t%s\nsearch_lnl\t%s\nmodel_file\t%s.best_model.nex\nstart_tree\t%s.treefile\n' \
  "$best_name" "$best_lnl" "$best_prefix" "$best_prefix" > "$OUT/params-from.txt"

{
  echo "started: $(date -Is)"
  echo "host: $(hostname)"
  echo "slurm_job: ${SLURM_JOB_ID:-none}"
  echo "model: fixed at $best_prefix.best_model.nex (highest-likelihood search: $best_name)"
  "$IQTREE" --version | head -1
} > "$OUT/environment.txt.tmp"

start=$(date +%s)
"$IQTREE" -s "$ALN" -st DNA -p "$best_prefix.best_model.nex" -T "$THREADS" --seed "$SEED" --keep-ident \
  -t "$best_prefix.treefile" --tree-fix --trees "$OUT/candidates.nwk" \
  --test 10000 --test-au --test-weight --prefix "$OUT/au" > "$OUT/au.stdout" 2> "$OUT/au.stderr"
printf 'au\t%s\t%s\n' "$(( $(date +%s) - start ))" "${SLURM_JOB_ID:-none}" > "$OUT/stage-seconds.tsv"
echo "finished: $(date -Is)" >> "$OUT/environment.txt.tmp"
mv "$OUT/environment.txt.tmp" "$OUT/environment.txt"
