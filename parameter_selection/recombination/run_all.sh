#!/usr/bin/env bash
# Run the whole tutorial from an empty folder: about 3 minutes on a laptop.
# Requires zipstrain 2.0.0, bowtie2, samtools, wgsim, curl and matplotlib
# (see the "Recombination mosaics" page of the ZipStrain documentation).
set -euo pipefail
cd "$(dirname "$0")"
threads=${THREADS:-4}

echo "== 1. Download the S. pneumoniae TIGR4 reference (NC_003028.3)"
mkdir -p data/reference
if [ ! -s data/reference/TIGR4.fna ]; then
  curl -sSfL -o data/reference/TIGR4.fna \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_003028.3&rettype=fasta&retmode=text"
fi

echo "== 2. Build strain B with planted mosaic blocks"
python scripts/01_make_strains.py --reference data/reference/TIGR4.fna --outdir data

echo "== 3. Simulate 30x reads for samples A, B and a 70/30 mix"
bash scripts/02_simulate_reads.sh data

echo "== 4. Map reads to strain A"
zipstrain map -i data/reads.csv -o mapped -f data/strain_A.fna -s data/reference.stb -t "$threads"

echo "== 5. Profile"
zipstrain profile -i mapped/samples.txt -f data/strain_A.fna -s data/reference.stb -r profiles --threads "$threads"

echo "== 6. Find mosaic regions and identical stretches"
python scripts/03_find_mosaics.py --profiles profiles --truth data/truth/planted_blocks.tsv --outdir results

echo "== 7. Draw the figures"
python scripts/04_plot_figures.py --truth data/truth/planted_blocks.tsv --results results --profiles profiles --outdir figures

echo "== 8. Export the A-vs-B mismatches for interactive/mosaic_explorer.html"
python scripts/05_export_explorer_data.py --profiles profiles --truth data/truth/planted_blocks.tsv \
  --reference data/strain_A.fna --out interactive/tutorial_data.js
