#!/usr/bin/env bash
# This shell script is used to correctly transfer the hg19 format of the granja data to the hg38 format. To be used before running any processing.

set -euo pipefail

in_dir="."
out_dir="processed"
mkdir -p "$out_dir"
chain="hg19ToHg38.over.chain.gz"

for f in "$in_dir"/granja_*_CD34_*.fragments.tsv.gz; do
  name=$(basename "$f" .fragments.tsv.gz)
  echo "=== $name ==="

  tmp=$(mktemp --suffix=.bed)

  echo "  [1/3] Decompress → BED6"
  zcat "$f" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' > "$tmp"

  echo "  [2/3] liftOver"
  liftOver -bedPlus=3 "$tmp" "$chain" \
            "$tmp.lift"  "$out_dir/${name}_unmapped_hg19.bed"

  echo "  [2b/3] Sorting for tabix"
  sort -k1,1 -k2,2n "$tmp.lift" -o "$tmp.lift.sorted"

  echo "  [3/3] Compressing & indexing"
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' "$tmp.lift.sorted" \
      | bgzip > "$out_dir/${name}_hg38.fragments.tsv.gz"
  tabix -p bed -0 -s 1 -b 2 -e 3 "$out_dir/${name}_hg38.fragments.tsv.gz"

  rm "$tmp" "$tmp.lift" "$tmp.lift.sorted"
  echo "✓ Done $name"
done

echo "All liftOvers complete – output in $out_dir/"