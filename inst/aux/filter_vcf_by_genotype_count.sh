#!/bin/bash

# Usage:
#   cat input.vcf | ./filter_vcf_by_genotype_count.sh <min_count> <min_classes> > output.vcf
#   ./filter_vcf_by_genotype_count.sh <min_count> <min_classes> < input.vcf > output.vcf
#
# Arguments:
#   min_count   - Minimum number of samples per genotype class (integer)
#   min_classes - Minimum number of genotype classes (2 or 3) that must meet min_count

MIN_COUNT="$1"
MIN_CLASSES="${2:-3}"

# Check inputs
if [[ -z "$MIN_COUNT" || -z "$MIN_CLASSES" || "$MIN_COUNT" =~ [^0-9] || "$MIN_CLASSES" =~ [^0-9] || "$MIN_CLASSES" -lt 2 || "$MIN_CLASSES" -gt 3 ]]; then
  echo "Usage: $0 <min_genotype_count> <min_genotype_classes>" >&2
  echo "  <min_genotype_count>: Minimum samples per genotype class (integer)" >&2
  echo "  <min_genotype_classes>: Number of genotype classes (2 or 3)" >&2
  echo "Example: cat input.vcf | $0 3 3 > filtered.vcf" >&2
  exit 1
fi

awk -v min="$MIN_COUNT" -v classes="$MIN_CLASSES" '
BEGIN {
  FS = OFS = "\t"
  kept = removed = 0
}
/^##/ {
  print; next
}
/^#CHROM/ {
  print; next
}
{
  g00 = g01 = g11 = 0

  for (i = 10; i <= NF; i++) {
    gt = $i
    if (gt == "." || gt == "./." || gt == ".|.") continue
    gsub("\\|", "/", gt)
    if (gt == "0/0") g00++
    else if (gt == "0/1" || gt == "1/0") g01++
    else if (gt == "1/1") g11++
  }

  class_pass = 0
  if (g00 >= min) class_pass++
  if (g01 >= min) class_pass++
  if (g11 >= min) class_pass++

  if (class_pass >= classes) {
    print
    kept++
  } else {
    removed++
  }
}
END {
  print "### Summary: " kept " variants kept, " removed " removed (threshold = " min ", classes = " classes ")" >> "/dev/stderr"
}
'