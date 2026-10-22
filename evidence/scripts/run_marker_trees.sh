#!/usr/bin/env bash
# The BAC120 marker trees for a cohort: the primary, its model baseline and its program baseline.
#
# The primary is FastTree 2.2.0 WAG+gamma through `gtdbtk infer`, so the parameterisation is GTDB's
# rather than ours (SET-15), on the untrimmed occupancy-filtered alignment (SET-16).
#
# The other three exist because SET-18 rule 3 requires every conflict figure to be quoted with a
# within-method baseline measured in the same statistic, and the pilot's baselines cannot simply be
# carried over:
#
#   * **Model choice** -- LG and JTT on the identical alignment. One parameter moved, alignment,
#     program and algorithm held. This is the exact analogue of Verticall's distance mode and
#     Mashtree's seed, and it is what makes the marker arm comparable to the other two.
#   * **Program choice** -- VeryFastTree 4.0.5 on the identical alignment, invoked with the flags
#     GTDB-Tk itself passes (`-wag -gamma`, support left on) plus `-double-precision`, which is
#     required to match the GTDB-Tk FastTree build: on an alignment this full of near-zero branches
#     single precision is not a cosmetic difference. The pilot's program-choice baseline was
#     FastTree against IQ-TREE; SET-15 removed IQ-TREE from the method, and a 1,090-tip partitioned
#     ML run is not warranted merely to produce a baseline.
#
# WAG stays primary for GTDB comparability even though Task 07 ranked it the weakest of FastTree's
# three matrices on this alignment: all three gave identical DEC-004 clades and no conflicting
# supported clade, so the choice costs nothing and buys method-matching.
#
# Usage: scripts/run_marker_trees.sh --run-dir DIR --alignment MSA [--threads N]
#                                    [--models "WAG LG JTT"] [--skip-veryfasttree]
set -euo pipefail

PROJECT=${PROJECT_ROOT}
RUN_DIR=""
MSA=""
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MODELS="WAG LG JTT"
SKIP_VFT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir)           RUN_DIR="$2"; shift 2 ;;
    --alignment)         MSA="$2"; shift 2 ;;
    --threads)           THREADS="$2"; shift 2 ;;
    --models)            MODELS="$2"; shift 2 ;;
    --skip-veryfasttree) SKIP_VFT=1; shift ;;
    *) echo "Error: unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -n "$RUN_DIR" ]] || { echo "Error: --run-dir is required" >&2; exit 2; }
[[ -s "$MSA" ]]     || { echo "Error: --alignment must name a non-empty file" >&2; exit 2; }

cd "$PROJECT"
module load gtdb-tk/2.7.2+232
export GTDBTK_DATA_PATH=/scratch/opendata/genomics/GTDB/releases/release232/
export PATH="/sw/local/rocky9/noarch/rcc/software/gtdb-tk/gtdbtk_2.7.2/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS="$THREADS"

LOGS="$RUN_DIR/logs"
mkdir -p "$LOGS" "$RUN_DIR/fasttree" "$RUN_DIR/veryfasttree"
MSA_SHA="$(sha256sum "$MSA" | cut -d' ' -f1)"

{
  echo "alignment: $MSA"
  echo "alignment_sha256: $MSA_SHA"
  echo "fasttree: $(FastTree 2>&1 | head -1)"
  echo "veryfasttree: $($PROJECT/.work/envs/veryfasttree/bin/VeryFastTree 2>&1 | sed -n 2p)"
  echo "gtdbtk: $(gtdbtk --version 2>&1 | head -1)"
  echo "threads: $THREADS"
  echo "started: $(date -Is)"
} > "$LOGS/marker-trees-environment.txt"

for model in $MODELS; do
  label="$(echo "$model" | tr '[:upper:]' '[:lower:]')"
  [[ "$model" == "WAG" ]] && label="filtered" || label="filtered-$label"
  out="$RUN_DIR/fasttree/$label"
  if [[ -s "$out/gtdbtk.unrooted.tree" ]]; then echo "$label already done"; continue; fi
  rm -rf "$out"; mkdir -p "$out"
  echo "== gtdbtk infer --prot_model $model --gamma -> $out"
  /usr/bin/time -v -o "$LOGS/fasttree-$label.time" gtdbtk infer \
    --msa_file "$MSA" --out_dir "$out" --prot_model "$model" --gamma --cpus "$THREADS" \
    > "$LOGS/fasttree-$label.stdout" 2> "$LOGS/fasttree-$label.stderr"
done

if [[ $SKIP_VFT -eq 0 ]]; then
  out="$RUN_DIR/veryfasttree/filtered"
  if [[ ! -s "$out/veryfasttree.nwk" ]]; then
    rm -rf "$out"; mkdir -p "$out"
    echo "== VeryFastTree -wag -gamma -double-precision -> $out"
    # The same flags GTDB-Tk passes FastTree (external/fasttree.py: -wag, -gamma, -log, support
    # left on), plus -double-precision to match the GTDB-Tk build's arithmetic.
    /usr/bin/time -v -o "$LOGS/veryfasttree.time" \
      "$PROJECT/.work/envs/veryfasttree/bin/VeryFastTree" \
        -wag -gamma -double-precision -threads "$THREADS" \
        -log "$out/veryfasttree.log" "$MSA" \
        > "$out/veryfasttree.nwk" 2> "$LOGS/veryfasttree.stderr"
  else
    echo "veryfasttree already done"
  fi
fi

echo "finished: $(date -Is)" >> "$LOGS/marker-trees-environment.txt"
find "$RUN_DIR/fasttree" "$RUN_DIR/veryfasttree" -name '*.tree' -o -name '*.nwk' \
  | sort | xargs sha256sum > "$RUN_DIR/marker-trees.sha256"
echo "marker trees complete"
