#!/usr/bin/env bash
# Exact genome size and G+C (mol%) of type whole assemblies from frozen local FASTA (read-only).
# Revised 9 Oct 2026 (coordinator review R1): integer A/C/G/T counts per contig from
# `seqkit fx2tab --base-count` (-C), summed over contigs before dividing. The previous
# version rebuilt counts from rounded per-contig percentages (-B) and is kept in
# superseded/2026-10-08-pre-review/scripts/.
# Usage: type_genome_gc.sh OUT.tsv ACC=PATH ...
# Per-contig counts are written to $(dirname OUT.tsv)/base-counts/ACC.base-counts.tsv.
set -euo pipefail
out=$1; shift
pc=$(dirname "$out")/base-counts
mkdir -p "$pc"
printf "accession\tfasta_path\tfasta_sha256\tn_sequences\ttotal_length\tA\tC\tG\tT\tgc_mol_percent_acgt\tseqkit_version\tnon_acgt_bases\tper_contig_counts\n" > "$out"
v=$(seqkit version | awk '{print $2}')
for kv in "$@"; do
  acc=${kv%%=*}; p=${kv#*=}
  sha=$(sha256sum "$p" | cut -d' ' -f1)
  f=$pc/$acc.base-counts.tsv
  { printf "name\tlength\tA\tC\tG\tT\n"; seqkit fx2tab -n -l -C A -C C -C G -C T "$p"; } > "$f"
  awk -F'\t' -v acc="$acc" -v p="$p" -v sha="$sha" -v v="$v" -v f="$f" '
    NR>1 {n++; L+=$2; a+=$3; c+=$4; g+=$5; t+=$6}
    END{printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%s\t%d\t%s\n", acc,p,sha,n,L,a,c,g,t,100*(g+c)/(a+c+g+t),v,L-(a+c+g+t),f}' "$f" >> "$out"
done
