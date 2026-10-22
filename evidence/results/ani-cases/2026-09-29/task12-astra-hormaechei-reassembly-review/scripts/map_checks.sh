#!/usr/bin/env bash
set -euo pipefail
export PATH=${USER_HOME}/sw/miniforge3/envs/autocycler/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=${PROJECT_ROOT}
R=$ROOT/results/ani-cases/2026-09-29/task12-astra-hormaechei-reassembly-review
D=$ROOT/data/reads/hormaechei-type-SAMN16357575
cd "$R"
{ minimap2 --version; samtools --version; seqkit version; python -c 'import pysam; print(pysam.__version__)'; } > logs/versions.txt
seqkit seq --min-len 20000 --threads 1 "$D/ccs/SRR15075454.hifi_reads.fastq" > work/hifi-ge20k.fq
seqkit seq --min-len 20000 --threads 1 "$D/SRR15075454.fastq.gz" > work/subreads-ge20k.fq
# -P retains all chains; map each model alone as a separate sensitivity below.
minimap2 -t 4 -cx map-hifi -P --cs=long work/junction-models.fa work/hifi-ge20k.fq > work/span-hifi-allchains.paf 2> logs/span-hifi.log
minimap2 -t 4 -cx map-pb -P --cs=long work/junction-models.fa work/subreads-ge20k.fq > work/span-subreads-allchains.paf 2> logs/span-subreads.log
for a in final published collapsed; do
 samtools faidx "work/$a.fa"
 minimap2 -t 4 -ax map-hifi --secondary=no "work/$a.fa" "$D/ccs/SRR15075454.hifi_reads.fastq" 2> "logs/map-$a.log" | samtools sort -@ 1 -m 1G -o "work/$a.bam" -
 samtools index "work/$a.bam"
done
# Preserve independent reproduction of the original counter without editing its script.
for a in final published; do
 samtools mpileup -B -Q 0 -q 0 -d 0 -f "work/$a.fa" "work/$a.bam" 2> "logs/pileup-$a.log" | python "$ROOT/results/ani-cases/2026-09-29/task12-opus55-hormaechei-reassembly/pileup_checks.py" --label "$a" > "outputs/opus-counter-$a.tsv"
done
