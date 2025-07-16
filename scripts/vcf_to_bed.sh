bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/ANN[\t%GT]\n' merged_severus.simpleann.vcf | \
awk -v OFS='\t' '
BEGIN {
  # Short sample names in order
  split("2436A_run3 2436A_run5 3188 6265D 6266D 6269C 6277A 6285B 6298B 6314E 6324B 6365A 6414C 6423A 6436B 6439C 6478A 6483C 6500D 6605D", sample_short);
}
{
  chrom = $1;
  start = $2 - 1;
  end = $3;
  svtype = $4;

  # Get first ANN entry and extract the functional annotation
  split($5, ann_parts, ",");
  split(ann_parts[1], ann_fields, "|");
  func_annot = ann_fields[2];  # e.g., exon_loss_variant

  for (i = 6; i <= NF; i++) {
    sample = sample_short[i-5];
    gt = $i;
    if (gt != "./.") {
      print chrom, start, end, sample, svtype, gt, func_annot;
    }
  }
}
' > simple-sv_annot.bed

