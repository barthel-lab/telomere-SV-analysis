#!/bin/bash
#SBATCH --job-name=snpeff_merged
#SBATCH --output=snpeff_merged_%j.out
#SBATCH --error=snpeff_merged_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00

module load Java/11.0.2

java -Xmx16g -jar ~/softwares/snpEff/snpEff.jar \
  -c ~/softwares/snpEff/snpEff.config -v CHM13v2.0 \
  /scratch/mjehangir/Glioma_project/severus/severus_out/sv_annotation/merged_severus_cohort_rename.vcf.gz > \
  /scratch/mjehangir/Glioma_project/severus/severus_out/sv_annotation/merged_severus_cohort.snpeff_annotation.vcf

