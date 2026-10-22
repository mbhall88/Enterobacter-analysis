#!/usr/bin/env bash
# Task 17 lines 4 and 5: gene trees and concordance factors on the named branches of one panel.
#
# IQ-TREE 3.1.1 is a testing instrument here; SET-15's FastTree tree is not replaced.
#
#   1. one gene tree per locus of --loci (IQ-TREE -S), ModelFinder per locus over a stated set;
#   2. for every reference tree in --references (name=path, repeatable): gene concordance
#      (--gcf) and, separately, likelihood site concordance (--scfl) -- IQ-TREE 3 refuses the two
#      in one call -- with the reference FIXED (--tree-fix). Task 07 recorded what happens without
#      it: IQ-TREE treats -t as a starting tree, searches, and annotates a different topology.
#   3. sCF is repeated with several seeds so its Monte-Carlo stability is measured, not assumed.
#
# Reference trees are the unconstrained ML tree and the constrained ML trees from
# run_topology_tests.sh on the same panel, so every competing resolution's branch exists in at
# least one reference.
#
# Usage: bash scripts/run_concordance.sh --alignment F --loci F --out DIR --model M
#          --references name=tree [--references ...] [--seqtype AA|DNA] [--threads N]
#          [--quartets N] [--seeds "1 2 3"]
set -euo pipefail

PROJECT=${PROJECT_ROOT}
cd "$PROJECT"
IQTREE=${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

ALN="" LOCI="" OUT="" MODEL="" PART="" SEQTYPE="AA" QUARTETS=1000 SEEDS="20260923 20260924 20260925"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
REFS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --alignment) ALN="$2"; shift 2 ;;
    --loci) LOCI="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --partition) PART="$2"; shift 2 ;;
    --references) REFS+=("$2"); shift 2 ;;
    --seqtype) SEQTYPE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --quartets) QUARTETS="$2"; shift 2 ;;
    --seeds) SEEDS="$2"; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -s "$ALN" && -s "$LOCI" && -n "$OUT" && ( -n "$MODEL" || -n "$PART" ) && ${#REFS[@]} -gt 0 ]] \
  || { echo "need --alignment, --loci, --out, --model or --partition, and --references" >&2; exit 2; }
# sCF is computed under the same model as the reference trees: a single -m for the markers, the
# codon-position partition (-p) for core genes.
if [[ -n "$PART" ]]; then scf_model=(-p "$PART"); else scf_model=(-m "$MODEL"); fi
if [[ -e "$OUT" ]]; then echo "refusing to overwrite $OUT" >&2; exit 2; fi
mkdir -p "$OUT/logs"

run() {
  local label="$1"; shift
  local start=$(date +%s)
  /usr/bin/time -v -o "$OUT/logs/$label.time" "$IQTREE" "$@" \
    > "$OUT/logs/$label.stdout" 2> "$OUT/logs/$label.stderr"
  echo -e "$label\t$(( $(date +%s) - start ))" >> "$OUT/stage-seconds.tsv"
}

{
  echo "iqtree: $("$IQTREE" --version | head -1)"
  echo "host: $(hostname)"; echo "slurm_job_id: ${SLURM_JOB_ID:-none}"; echo "threads: $THREADS"
  echo "alignment: $ALN $(sha256sum "$ALN" | cut -d' ' -f1)"
  echo "loci: $LOCI $(sha256sum "$LOCI" | cut -d' ' -f1) ($(grep -c 'charset' "$LOCI") loci)"
  echo "sCF model: ${MODEL:-partition $PART}"; echo "sCF quartets: $QUARTETS"; echo "sCF seeds: $SEEDS"
  for ref in "${REFS[@]}"; do echo "reference: $ref $(sha256sum "${ref#*=}" | cut -d' ' -f1)"; done
  echo "started: $(date -Is)"
} > "$OUT/environment.txt"

if [[ "$SEQTYPE" == AA ]]; then
  gene_models=(-m MFP --mset LG,WAG,JTT,Q.PFAM,Q.INSECT,Q.PLANT,Q.YEAST --mrate E,G,I+G)
else
  gene_models=(-m MFP --mset JC,K2P,HKY,TN,GTR --mrate E,G,I+G)
fi
run loci -s "$ALN" -st "$SEQTYPE" -S "$LOCI" "${gene_models[@]}" -T "$THREADS" \
  --seed 20260923 --keep-ident --prefix "$OUT/loci"

# Gene trees from near-identical sequences resolve zero-length branches arbitrarily, and gCF counts
# whatever was written. A second gCF is therefore computed on gene trees with every internal branch
# shorter than one expected substitution across that gene collapsed (SET-18 rule 10's
# one-unit-of-evidence rule, in each gene's own units). Both are reported; the collapsed one is
# the one read.
python3 - "$OUT/loci.treefile" "$LOCI" "$OUT/loci.collapsed.treefile" <<'PY'
import re, sys
sys.path.insert(0, "scripts")
from summarise_concordance import collapse_short
trees = [l.strip() for l in open(sys.argv[1]) if l.strip()]
loci = [(int(a), int(b)) for a, b in re.findall(r"charset\s+\S+\s*=\s*(\d+)-(\d+)\s*;", open(sys.argv[2]).read())]
if len(trees) != len(loci):
    sys.exit(f"{len(trees)} gene trees for {len(loci)} loci")
with open(sys.argv[3], "w") as out:
    for tree, (start, end) in zip(trees, loci):
        out.write(collapse_short(tree, 1.0 / (end - start + 1)) + "\n")
PY

for ref in "${REFS[@]}"; do
  name="${ref%%=*}" tree="${ref#*=}"
  run "gcf-$name" -t "$tree" --gcf "$OUT/loci.treefile" --df-tree --cf-verbose \
    -T 1 --prefix "$OUT/$name.gcf"
  run "gcfc-$name" -t "$tree" --gcf "$OUT/loci.collapsed.treefile" --cf-verbose \
    -T 1 --prefix "$OUT/$name.gcfc"
  for seed in $SEEDS; do
    run "scf-$name-$seed" -t "$tree" --tree-fix -s "$ALN" -st "$SEQTYPE" "${scf_model[@]}" \
      --scfl "$QUARTETS" --seed "$seed" -T "$THREADS" --keep-ident --prefix "$OUT/$name.scf-$seed"
    if grep -q "BETTER TREE FOUND" "$OUT/$name.scf-$seed.log"; then
      echo "tree search ran despite --tree-fix for $name" >&2; exit 4
    fi
  done
done

echo "finished: $(date -Is)" >> "$OUT/environment.txt"
echo "concordance complete: $OUT"
