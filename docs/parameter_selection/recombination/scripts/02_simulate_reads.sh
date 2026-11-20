#!/usr/bin/env bash
# Simulate 30x paired-end Illumina-like reads for three samples:
#   sample_A  reads from strain A
#   sample_B  reads from strain B
#   sample_M  a mixed population: 70% strain B + 30% strain A
# Usage: scripts/02_simulate_reads.sh <data-dir>
set -euo pipefail

data=${1:-data}
reads="$data/reads"
mkdir -p "$reads"

genome_length=$(grep -v '>' "$data/strain_A.fna" | tr -d '\n' | wc -c | tr -d ' ')
read_length=150
coverage=30
pairs=$(( genome_length * coverage / (2 * read_length) ))

# -e base error rate, -r 0 / -R 0 / -X 0: wgsim adds no mutations of its own,
# so every difference between samples is one we planted.
simulate() {  # simulate <genome> <pairs> <seed> <prefix>
  wgsim -N "$2" -1 "$read_length" -2 "$read_length" -d 400 -s 50 \
        -e 0.001 -r 0 -R 0 -X 0 -S "$3" "$1" "$4_1.fq" "$4_2.fq" > /dev/null
}

simulate "$data/strain_A.fna" "$pairs" 11 "$reads/sample_A"
simulate "$data/strain_B.fna" "$pairs" 22 "$reads/sample_B"

# Mixed population: same total depth, 70% B and 30% A.
simulate "$data/strain_B.fna" $(( pairs * 7 / 10 )) 33 "$reads/mixB"
simulate "$data/strain_A.fna" $(( pairs * 3 / 10 )) 44 "$reads/mixA"
for mate in 1 2; do
  cat "$reads/mixB_$mate.fq" "$reads/mixA_$mate.fq" > "$reads/sample_M_$mate.fq"
  rm "$reads/mixB_$mate.fq" "$reads/mixA_$mate.fq"
done

gzip -f "$reads"/*.fq

{
  echo "sample_name,reads1,reads2"
  for sample in sample_A sample_B sample_M; do
    echo "$sample,$reads/${sample}_1.fq.gz,$reads/${sample}_2.fq.gz"
  done
} > "$data/reads.csv"

echo "Simulated $pairs read pairs (${coverage}x) per sample into $reads"
