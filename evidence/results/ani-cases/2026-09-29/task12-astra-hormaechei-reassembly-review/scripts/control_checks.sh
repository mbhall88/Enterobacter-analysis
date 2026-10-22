#!/usr/bin/env bash
set -euo pipefail
set -x
R=results/ani-cases/2026-09-29/task12-astra-hormaechei-reassembly-review
export PATH=${USER_HOME}/sw/miniforge3/envs/autocycler/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$R/scripts/fixtures.py"
minimap2 -t 2 -cx map-pb -P --cs=long "$R/work/junction-models.fa" "$R/work/span-controls.fa" > "$R/work/span-controls.paf" 2> "$R/logs/span-controls.log"
for a in final collapsed; do
  minimap2 -t 2 -ax map-hifi --secondary=no "$R/work/$a.fa" "$R/work/repeat-power.fq" 2> "$R/logs/$a-sim-map.log" | samtools sort -@ 0 -m 1G -o "$R/work/$a-sim.bam"
  samtools index "$R/work/$a-sim.bam"
done
python "$R/scripts/read_audit.py" --sim
for a in GCA_001875655.1 GCA_000213995.1; do
  minimap2 -t 2 -cx asm5 -P --cs=long "$R/work/junction-models.fa" "$R/type-census/raw/${a}_genomic.fna.gz" > "$R/work/$a-models.paf" 2> "$R/logs/$a-models.log"
done
