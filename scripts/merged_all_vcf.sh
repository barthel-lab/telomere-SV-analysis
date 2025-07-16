#!/bin/bash

parent_dir="/scratch/mjehangir/Glioma_project/severus/severus_out"
merged_vcf="${parent_dir}/merged_severus_cohort.vcf.gz"

# Step 1: Find all Severus VCF files
vcf_files=$(find "$parent_dir" -type f -name "severus_somatic.vcf")

if [ -z "$vcf_files" ]; then
    echo "No Severus VCF files found!"
    exit 1
fi

echo "Found Severus VCF files:"
echo "$vcf_files"

# Step 2: Compress and index each VCF if not already done
while IFS= read -r vcf; do
    if [[ ! -f "${vcf}.gz" ]]; then
        echo "Compressing $vcf"
        bgzip -c "$vcf" > "${vcf}.gz"
    fi
    if [[ ! -f "${vcf}.gz.tbi" ]]; then
        echo "Indexing ${vcf}.gz"
        tabix -p vcf "${vcf}.gz"
    fi
done <<< "$vcf_files"

# Step 3: Prepare list of compressed VCFs for merging
compressed_vcfs=""
while IFS= read -r vcf; do
    compressed_vcfs+=" ${vcf}.gz"
done <<< "$vcf_files"

# Step 4: Merge VCFs into one multi-sample VCF
echo "Merging all VCFs into $merged_vcf"
bcftools merge $compressed_vcfs -Oz -o "$merged_vcf"

# Step 5: Index merged VCF
tabix -p vcf "$merged_vcf"

echo "Merged VCF created at: $merged_vcf"

