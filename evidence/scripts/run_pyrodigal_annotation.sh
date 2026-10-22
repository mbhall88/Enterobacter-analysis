#!/usr/bin/env bash
# Task 17: call genes with pyrodigal on every genome of the focused subsets, for Panaroo.
#
# One gene call per genome, shared by all three conflicts, so a genome that sits in two subsets is
# annotated once and identically. Only gene structure matters here: orthologs are clustered by
# sequence, so no functional annotation (Bakta/Prokka) is needed and none is added, which keeps
# every gene call from one pinned program with no database behind it.
#
# Each derivative is re-hashed against the subset manifest before it is read.
# Outputs <tip>.gff, <tip>.faa, <tip>.ffn, and a list file of "<tip>.gff <derivative>" pairs in the
# form Panaroo 1.8.0 accepts (its __main__.py reads two whitespace-separated paths per line).
#
# Usage: bash scripts/run_pyrodigal_annotation.sh <inputs.tsv> <outdir>
set -euo pipefail

INPUTS="${1:?usage: run_pyrodigal_annotation.sh <inputs.tsv> <outdir>}"
OUT="${2:?missing outdir}"
PROJECT=${PROJECT_ROOT}
cd "$PROJECT"
ENV="$PROJECT/.work/envs/panaroo"
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-4}}"

if [[ -e "$OUT/gff" ]]; then
  echo "refusing to overwrite existing $OUT/gff" >&2; exit 2
fi
mkdir -p "$OUT/gff" "$OUT/faa" "$OUT/ffn" "$OUT/logs"
"$ENV/bin/pyrodigal" -V > "$OUT/pyrodigal-version.txt" 2>&1

annotate() {
  local tip="$1" fasta="$2" expected="$3" out="$4" env="$5"
  local got
  got=$(sha256sum "$fasta" | cut -d' ' -f1)
  if [[ "$got" != "$expected" ]]; then
    echo "$tip: $fasta hash $got does not match manifest $expected" >&2; return 3
  fi
  "$env/bin/pyrodigal" -i "$fasta" -f gff -o "$out/gff/$tip.gff" -a "$out/faa/$tip.faa" \
    -d "$out/ffn/$tip.ffn" -g 11 -p single > "$out/logs/$tip.log" 2>&1
  printf '%s\t%s\n' "$tip" "$(grep -c $'\tCDS\t' "$out/gff/$tip.gff")"
}
export -f annotate

start=$(date +%s)
tail -n +2 "$INPUTS" \
  | xargs -P "$CPUS" -L 1 bash -c 'annotate "$0" "$1" "$2" "'"$OUT"'" "'"$ENV"'"' \
  > "$OUT/cds-counts.unsorted.tsv"
sort "$OUT/cds-counts.unsorted.tsv" > "$OUT/cds-counts.tsv" && rm "$OUT/cds-counts.unsorted.tsv"
end=$(date +%s)

expected=$(( $(wc -l < "$INPUTS") - 1 ))
done_count=$(wc -l < "$OUT/cds-counts.tsv")
if [[ "$done_count" -ne "$expected" ]]; then
  echo "annotated $done_count of $expected genomes" >&2; exit 4
fi
tail -n +2 "$INPUTS" | while IFS=$'\t' read -r tip fasta _; do
  printf '%s %s\n' "$PROJECT/$OUT/gff/$tip.gff" "$PROJECT/$fasta"
done > "$OUT/panaroo-inputs-all.txt"
( cd "$OUT" && sha256sum gff/*.gff faa/*.faa > annotation.sha256 )
echo "annotated $done_count genomes in $((end - start)) s at $CPUS CPUs"
