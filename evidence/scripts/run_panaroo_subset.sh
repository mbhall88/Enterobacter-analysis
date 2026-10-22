#!/usr/bin/env bash
# Task 17 line 5: cluster orthologs and align the core genes of one frozen focused subset.
#
# Panaroo 1.8.0 (pinned at .work/envs/panaroo) on the pyrodigal gene calls from
# run_pyrodigal_annotation.sh. The subset's ids file decides the input set; every id must have
# exactly one GFF/FASTA pair, and the isolates Panaroo reports are checked against the ids file
# afterwards, so a silently dropped genome fails the run.
#
# Settings, all stated so none is a hidden default:
#   --clean-mode strict      Panaroo's recommended mode for phylogenetic core-genome work
#   --remove-invalid-genes   drop calls Panaroo cannot translate cleanly, counted in the log
#   -a core --aligner mafft --codons
#                            per-gene codon alignments built from protein alignments, so the
#                            genus-level divergence of C1/C2 (~88-92% ANI) aligns in frame
#   --core_threshold 0.95    genes in >=95% of the subset; per-gene analyses may filter further
#   -c 0.98 -f 0.7           Panaroo defaults, written out
#
# Usage: bash scripts/run_panaroo_subset.sh <subset.ids.txt> <panaroo-inputs-all.txt> <outdir>
set -euo pipefail

IDS="${1:?usage: run_panaroo_subset.sh <ids> <inputs-all> <outdir>}"
ALL="${2:?missing inputs list}"
OUT="${3:?missing outdir}"
PROJECT=${PROJECT_ROOT}
cd "$PROJECT"
ENV="$PROJECT/.work/envs/panaroo"
export PATH="$ENV/bin:$PATH"
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-4}}"

if [[ -e "$OUT" ]]; then
  echo "refusing to overwrite existing $OUT" >&2; exit 2
fi
mkdir -p "$OUT"
LIST="$OUT/panaroo-inputs.txt"
while read -r tip; do
  [[ -n "$tip" ]] || continue
  line=$(grep -F "/gff/$tip.gff " "$ALL" || true)
  if [[ $(printf '%s\n' "$line" | grep -c .) -ne 1 ]]; then
    echo "$tip: expected exactly one input pair in $ALL" >&2; exit 3
  fi
  printf '%s\n' "$line"
done < "$IDS" > "$LIST"
[[ $(wc -l < "$LIST") -eq $(grep -c . "$IDS") ]] || { echo "input list incomplete" >&2; exit 3; }

{
  echo "host: $(hostname)"; echo "slurm_job_id: ${SLURM_JOB_ID:-none}"; echo "cpus: $CPUS"
  echo "panaroo: $(panaroo --version 2>&1)"; echo "mafft: $(mafft --version 2>&1 | head -1)"
  echo "cd-hit: $(cd-hit -h 2>&1 | head -1)"
  echo "ids: $IDS ($(grep -c . "$IDS") genomes, sha256 $(sha256sum "$IDS" | cut -d' ' -f1))"
  echo "started: $(date -Is)"
} > "$OUT/environment.txt"

/usr/bin/time -v -o "$OUT/panaroo.time" panaroo -i "$LIST" -o "$OUT/panaroo" \
  --clean-mode strict --remove-invalid-genes -c 0.98 -f 0.7 \
  -a core --aligner mafft --codons --core_threshold 0.95 -t "$CPUS" \
  > "$OUT/panaroo.stdout" 2> "$OUT/panaroo.stderr"

echo "finished: $(date -Is)" >> "$OUT/environment.txt"
python3 - "$IDS" "$OUT/panaroo/gene_presence_absence.csv" <<'EOF'
import csv, sys
ids = {l.strip() for l in open(sys.argv[1]) if l.strip()}
header = next(csv.reader(open(sys.argv[2])))
isolates = set(header[3:])
if isolates != ids:
    sys.exit(f"panaroo isolates differ from ids: missing {sorted(ids - isolates)[:5]}, "
             f"unexpected {sorted(isolates - ids)[:5]}")
print(f"panaroo isolates verified: {len(isolates)} == ids file")
EOF
echo "panaroo complete for $IDS"
