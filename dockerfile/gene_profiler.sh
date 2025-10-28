#!/bin/bash
set -euo pipefail

gene_table=$1
bam=$2
output_pileup=$3
threads=$4

unique_stamp=$(date +%s%N)_$$
mkdir -p tmp
tmpdir=$(mktemp -d "tmp/gene_pileup_$(basename "$bam" .bam)_${unique_stamp}_XXXXXX")

cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT
awk '{print $2":"$3"-"$4, $1}' "$gene_table" | \
parallel -j$4 --tmpdir "$tmpdir" --colsep ' ' --no-run-if-empty --halt now,fail=1 '
  region={1}
  gene={2}
  samtools mpileup -aa -A -r "$region" '"$bam"' | \
  awk -v gene="$gene" -v OFS="\t" '\''{
    gsub(/\^./,"",$5); gsub(/\$/,"",$5); gsub(/[+-][0-9]+[ACGTNacgtn]+/,"",$5);
    a=g=c=t=0;
    for(i=1;i<=length($5);i++){
      b=toupper(substr($5,i,1));
      if(b=="A") a++; else if(b=="C") c++; else if(b=="G") g++; else if(b=="T") t++;
    }
    print gene, $1, $2, a, c, g, t;
  }'\'' > '"$tmpdir"'/pileup_"$gene".tsv
'

# Combine all per-gene files into the final output
# Optional: add header
printf "gene\tchrom\tpos\tA\tC\tG\tT\n" > "$output_pileup"
cat "$tmpdir"/pileup_*.tsv >> "$output_pileup"
