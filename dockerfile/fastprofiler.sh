#!/bin/bash
set -euo pipefail

bedfile=$1
bam=$2
output_pileup=$3
gene_range_table=$4
threads=$5

unique_stamp=$(date +%s%N)_$$
mkdir -p tmp
tmpdir=$(mktemp -d "tmp/gene_pileup_$(basename "$bam" .bam)_${unique_stamp}_XXXXXX")

cleanup() {
    rm -rf "$tmpdir"
}


trap cleanup EXIT

nl=$(cat "$bedfile" | wc -l)
if [ "$nl" -lt $threads ]; then
    threads=$nl
fi
split -n l/$threads  "$bedfile" "$tmpdir"/scaffold_table_part_


ls "$tmpdir"/scaffold_table_part_* | \
parallel -j$threads --tmpdir "$tmpdir" --colsep ' ' --no-run-if-empty --halt now,fail=1 '
  batch={1}
  samtools mpileup -A -l "$batch" '"$bam"' | fast_profile utilities process_mpileup --gene-range-table-loc '"$gene_range_table"' --batch-bed "$batch" --output-file "$batch".parquet
'
fast_profile utilities merge_parquet --input-dir "$tmpdir"  --output-file "$output_pileup".parquet