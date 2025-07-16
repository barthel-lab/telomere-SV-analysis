# Define your short sample names in order (adjust as needed)
samples_short=("2436A_run3" "2436A_run5" "3188" "6265D" "6266D" "6269C" "6277A" "6285B" "6298B" "6314E" "6324B" "6365A" "6414C" "6423A" "6436B" "6439C" "6478A" "6483C" "6500D" "6605D")

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/ANN\t%INFO/SIMPLE_ANN[\t%GT]\n' merged_severus.simpleann.vcf | \
awk -v OFS='\t' -v samples="${samples_short[*]}" '
BEGIN {
  split(samples, sample_names, " ");
  print "CHROM", "POS", "END", "SVTYPE", "SnpEff_Effect", "SIMPLE_Effect", "SIMPLE_Label", "SIMPLE_Gene", "Sample", "Genotype";
}
{
  chrom=$1; pos=$2; end=$3; svtype=$4;

  # Parse SnpEff ANN for first effect
  split($5, ann_entries, ",");
  split(ann_entries[1], ann_fields, "|");
  split(ann_fields[2], effects, "&");
  snpeff_effect = effects[1];

  # Parse SIMPLE_ANN entries
  n=split($6, anns, ",");
  for(i=1; i<=n; i++) {
    split(anns[i], fields, "|");
    simple_effect = (length(fields) >= 2) ? fields[2] : "-";
    simple_gene   = (length(fields) >= 3) ? fields[3] : "-";
    simple_label  = (length(fields) >= 5) ? fields[5] : "-";

    # For each sample column (starting from $7)
    for(j=7; j<=NF; j++) {
      gt = $j;
      if(gt != "./." && gt != "." && gt != "0/0") {
        print chrom, pos, end, svtype, snpeff_effect, simple_effect, simple_label, simple_gene, sample_names[j-6], gt;
      }
    }
  }
}
' > sv_full_with_samples.tsv

