#!/bin/bash

# Set your main BAM directory
BAM_ROOT="/tgen_labs/barthel/fastq_DoNotTouch/Robert_Jenkins_Mayo/Dorado_bams"

# Output directory
OUTPUT_DIR="extracted_gene_bams"
mkdir -p "$OUTPUT_DIR"

# Loop through all BAMs
find "$BAM_ROOT" -name "*.sorted.bam" | while read bam; do
  sample=$(basename "$(dirname "$bam")")  # Gets folder name like 3188_merged_mod
  echo "Processing $sample..."

  while read chr start end gene; do
    outbam="${OUTPUT_DIR}/${sample}_${gene}.bam"
    samtools view -b "$bam" "${chr}:${start}-${end}" > "$outbam"
    samtools index "$outbam"
  done < merged_gene_regions.bed
done

