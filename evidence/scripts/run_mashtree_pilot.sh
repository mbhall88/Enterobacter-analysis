#!/usr/bin/env bash
# Task 09 -- matched Mashtree baseline on the frozen 388-tip panel.
#
# Reads input paths from a manifest into an argument array. It never globs
# MoniqueC_ECC_Study/all_assemblies/Enterobacter_Chr_Clean/: that directory holds a clinical query
# FASTA that is not a reference genome.
#
# Every sketch parameter is passed explicitly even where it equals the Mashtree default, so the
# command line is the provenance record and a later default change cannot silently alter the tree.
# Mashtree's sketching is seeded (mash sketch -S), so the seed is part of reproducibility.
#
# Usage:
#   scripts/run_mashtree_pilot.sh <manifest.list> <output-dir> [numcpus]
#
# Requires MASHTREE_ENV to point at the pinned prefix environment.
set -euo pipefail

manifest="${1:?usage: run_mashtree_pilot.sh <manifest.list> <output-dir> [numcpus]}"
outdir="${2:?usage: run_mashtree_pilot.sh <manifest.list> <output-dir> [numcpus]}"
numcpus="${3:-8}"

MASHTREE_ENV="${MASHTREE_ENV:-${PROJECT_ROOT}/.work/envs/mashtree}"
export PATH="${MASHTREE_ENV}/bin:${PATH}"

# Mashtree defaults, restated explicitly. --mindepth and --genomesize are read only for FASTQ
# input; every input here is an assembly, so they do not reach `mash sketch`. They are still
# recorded because "the default" is not a reproducible statement on its own.
KMERLENGTH=21
# Overridable only for the declared sketch-size sensitivity; the default is Mashtree's own and is
# what Task 09 and every seed replicate use, so leaving it unset reproduces those runs exactly.
SKETCH_SIZE="${MASHTREE_SKETCH_SIZE:-10000}"
SEED="${MASHTREE_SEED:-42}"   # overridable only for the declared seed-sensitivity run
GENOMESIZE=5000000
MINDEPTH=5
SORT_ORDER=ABC
TRUNCLENGTH=250
SIGFIGS=10

mkdir -p "${outdir}"
if [[ -e "${outdir}/mashtree.nwk" ]]; then
  echo "refusing to overwrite ${outdir}/mashtree.nwk" >&2
  exit 1
fi

# One derivative FASTA path per line -> argument array. No wildcard anywhere.
mapfile -t pilot_fastas < "${manifest}"
echo "inputs: ${#pilot_fastas[@]} from ${manifest}" >&2
for f in "${pilot_fastas[@]}"; do
  [[ -s "${f}" ]] || { echo "missing or empty input: ${f}" >&2; exit 1; }
done

# --tempdir keeps the sketches and the mash distance table for provenance instead of deleting
# them. --save-sketches is given an explicit directory: called bare it writes .msh files beside
# the source FASTAs, which here are Task 06's outputs and must not be written into.
{
  printf 'mashtree_version: %s\n' "$(mashtree --version 2>&1)"
  printf 'mash_version: %s\n' "$(mash --version 2>&1)"
  printf 'mash_path: %s\n' "$(command -v mash)"
  printf 'quicktree_version: %s\n' "$(quicktree -v 2>&1 | head -1)"
  printf 'numcpus: %s\n' "${numcpus}"
  printf 'kmerlength: %s\nsketch_size: %s\nseed: %s\n' "${KMERLENGTH}" "${SKETCH_SIZE}" "${SEED}"
  printf 'genomesize: %s\nmindepth: %s\n' "${GENOMESIZE}" "${MINDEPTH}"
  printf 'sort_order: %s\ntrunclength: %s\nsigfigs: %s\n' "${SORT_ORDER}" "${TRUNCLENGTH}" "${SIGFIGS}"
  printf 'inputs: %s\n' "${#pilot_fastas[@]}"
} > "${outdir}/run-parameters.txt"

/usr/bin/time -v mashtree \
  --numcpus "${numcpus}" \
  --kmerlength "${KMERLENGTH}" \
  --sketch-size "${SKETCH_SIZE}" \
  --seed "${SEED}" \
  --genomesize "${GENOMESIZE}" \
  --mindepth "${MINDEPTH}" \
  --sort-order "${SORT_ORDER}" \
  --truncLength "${TRUNCLENGTH}" \
  --sigfigs "${SIGFIGS}" \
  --tempdir "${outdir}/tmp" \
  --save-sketches "${outdir}/sketches" \
  --outmatrix "${outdir}/mash-distances.tsv" \
  --outtree "${outdir}/mashtree.nwk" \
  "${pilot_fastas[@]}" \
  > "${outdir}/mashtree.stdout" \
  2> "${outdir}/mashtree.time"

echo "wrote ${outdir}/mashtree.nwk" >&2
