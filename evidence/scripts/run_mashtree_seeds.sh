#!/usr/bin/env bash
# A recorded Mash seed sweep over any cohort, plus the majority-rule consensus SET-19 requires.
#
# Why a sweep at all. Task 09 changed one seed and 22.3% of splits moved; Task 10 ran twenty and
# found only 200 of 385 splits present in all of them. A single Mashtree draw is therefore not a
# usable comparator, and any Mashtree-versus-other-method figure has to be read against the seed
# distribution rather than against nothing.
#
# What the frequencies are. **Sketch stability** -- how reproducible a split is when mash re-draws
# its 10,000-hash sketch. SET-18 rules 2 and 6 forbid putting them in a column with, or on the
# scale of, the marker tree's bootstrap or SH-like support. They are not support.
#
# What SET-19 changed. The comparator is the majority-rule consensus of the replicates, not a
# filtered single draw: a filter applied to one draw can only remove splits, and 642 splits at
# pilot scale appeared in a replicate but not in the seed-42 tree, 25 of them in a majority.
# The single-seed tree is still built and kept, so already-circulated figures stay traceable.
#
# Usage:
#   scripts/run_mashtree_seeds.sh --list FILE --out DIR [--seeds "1001 1002 ..."]
#                                 [--primary-seed 42] [--cpus N] [--expected-tips FILE]
set -euo pipefail

ROOT=${PROJECT_ROOT}
LIST=""
OUT=""
SEEDS="1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1017 1018 1019 1020"
PRIMARY_SEED=42
CPUS="${SLURM_CPUS_PER_TASK:-8}"
EXPECTED=""
SKETCH_SIZE=10000
BUILD_PRIMARY=1
TREES_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --list)           LIST="$2"; shift 2 ;;
    --out)            OUT="$2"; shift 2 ;;
    --seeds)          SEEDS="$2"; shift 2 ;;
    --primary-seed)   PRIMARY_SEED="$2"; shift 2 ;;
    --cpus)           CPUS="$2"; shift 2 ;;
    --expected-tips)  EXPECTED="$2"; shift 2 ;;
    --sketch-size)    SKETCH_SIZE="$2"; shift 2 ;;
    --no-primary)     BUILD_PRIMARY=0; shift ;;
    --trees-only)     TREES_ONLY=1; shift ;;
    *) echo "Error: unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -n "$LIST" && -s "$LIST" ]] || { echo "Error: --list must name a non-empty file" >&2; exit 2; }
[[ -n "$OUT" ]] || { echo "Error: --out is required" >&2; exit 2; }

cd "$ROOT"
export MASHTREE_ENV="${MASHTREE_ENV:-$ROOT/.work/envs/mashtree}"
export MASHTREE_SKETCH_SIZE="$SKETCH_SIZE"
mkdir -p "$OUT/trees" "$OUT/logs"

# The single-seed draw, kept so already-circulated figures stay traceable (SET-19). The
# sketch-size sensitivity does not need it: it is compared consensus against consensus.
if [[ $BUILD_PRIMARY -eq 1 && ! -s "$OUT/trees/primary-seed-$PRIMARY_SEED.nwk" ]]; then
  dir="$OUT/work/primary-seed-$PRIMARY_SEED"
  rm -rf "$dir"
  MASHTREE_SEED="$PRIMARY_SEED" scripts/run_mashtree_pilot.sh "$LIST" "$dir" "$CPUS"
  cp "$dir/mashtree.nwk"       "$OUT/trees/primary-seed-$PRIMARY_SEED.nwk"
  cp "$dir/run-parameters.txt" "$OUT/logs/primary-seed-$PRIMARY_SEED.run-parameters.txt"
  grep -E 'Elapsed|Maximum resident|Percent of CPU|Exit status' "$dir/mashtree.time" \
    > "$OUT/logs/primary-seed-$PRIMARY_SEED.resources.txt"
  cp "$dir/mash-distances.tsv" "$OUT/primary-mash-distances.tsv"
  rm -rf "$dir"
fi

for seed in $SEEDS; do
  if [[ -s "$OUT/trees/seed-$seed.nwk" ]]; then echo "seed $seed already done"; continue; fi
  dir="$OUT/work/seed-$seed"
  rm -rf "$dir"
  MASHTREE_SEED="$seed" scripts/run_mashtree_pilot.sh "$LIST" "$dir" "$CPUS"
  cp "$dir/mashtree.nwk"       "$OUT/trees/seed-$seed.nwk"
  cp "$dir/run-parameters.txt" "$OUT/logs/seed-$seed.run-parameters.txt"
  grep -E 'Elapsed|Maximum resident|Percent of CPU|Exit status' "$dir/mashtree.time" \
    > "$OUT/logs/seed-$seed.resources.txt"
  # Sketches, tempdir and the NxN distance matrix are reproducible from the recorded seed and
  # parameters and are large; only the tree and its provenance are kept per replicate.
  rm -rf "$dir"
done

echo "seed distribution complete: $(ls "$OUT"/trees/seed-*.nwk | wc -l) replicate trees"

# --trees-only exists so a large sweep can be fanned across many array tasks writing into one
# output directory. Each task must NOT build a consensus of whatever subset happens to be present:
# that would race, and a consensus of a partial set is simply a different estimate. The consensus
# is built once, afterwards, over the complete set.
if [[ $TREES_ONLY -eq 1 ]]; then
  echo "--trees-only: stopping before the consensus"
  exit 0
fi

python3 scripts/majority_rule_consensus.py \
  $(for f in "$OUT"/trees/seed-*.nwk; do printf -- "--replicate %s " "$f"; done) \
  ${EXPECTED:+--expected-tips "$EXPECTED"} \
  --statistic-name mash_sketch_stability \
  --output-prefix "$OUT/consensus/majority-rule" \
  > "$OUT/logs/consensus.stdout" 2> "$OUT/logs/consensus.stderr"

# The seed-stable filter Task 10 used, rebuilt here only so the two can be compared directly.
if [[ $BUILD_PRIMARY -eq 1 ]]; then
python3 scripts/summarize_seed_stability.py \
  --reference "$OUT/trees/primary-seed-$PRIMARY_SEED.nwk" \
  $(for f in "$OUT"/trees/seed-*.nwk; do printf -- "--replicate %s " "$f"; done) \
  --threshold 0.95 --threshold 1.0 \
  --output-prefix "$OUT/seed-stability/seed-stability" \
  > "$OUT/logs/seed-stability.stdout" 2> "$OUT/logs/seed-stability.stderr"
fi

echo "consensus written under $OUT (sketch size $SKETCH_SIZE)"
