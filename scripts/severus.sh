#!/bin/bash
#SBATCH --job-name=sv
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

source ~/.bashrc

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "SAMPLE: $SAMPLE"

BASE_DIR="/tgen_labs/barthel/fastq_DoNotTouch/Robert_Jenkins_Mayo/Dorado_bams"
BAM="${BASE_DIR}/${SAMPLE}_merged_mod/${SAMPLE}_merged_mod.sorted.bam"

echo "BAM path: $BAM"

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found at $BAM" >&2
    exit 1
fi

PON="/scratch/mjehangir/Glioma_project/savana/Severus/pon/PoN_1000G_chm13.tsv.gz"
VNTR_BED="/scratch/mjehangir/Glioma_project/savana/Severus/vntrs/chm13.bed"
PHASING_VCF=""  # Provide path if available

mkdir -p severus_out/${SAMPLE} logs results

mamba activate severus_env

PHASING_OPTION=""
if [[ -n "$PHASING_VCF" ]]; then
    PHASING_OPTION="--phasing-vcf $PHASING_VCF"
fi

severus --target-bam "$BAM" \
        --out-dir "severus_out/${SAMPLE}" \
        -t 16 \
        $PHASING_OPTION \
        --vntr-bed "$VNTR_BED" \
        --PON "$PON" \
        > logs/${SAMPLE}_severus.out 2> logs/${SAMPLE}_severus.err

conda deactivate

