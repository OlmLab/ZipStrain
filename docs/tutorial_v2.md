# Tutorial

ZipStrain is designed to be flexible and can be run in a wide variety of ways. The three main decision points are:

| Decision | Option A | Option B |
|----------|----------|----------|
| **Execution** | Python command line interface (CLI) | Nextflow pipeline |
| **Reads** | Local reads (self-sequenced) | Public reads (SRA) |
| **Genome database** | Local genomes | Auto-built from public databases |

These three decisions can be mixed and matched in any combination. The tutorials below cover the most common use cases.

---

## Tutorial #1 — Python CLI with local reads

The following tutorial goes through an example run of ZipStrain using the Python CLI. You can follow along with your own data, or use a small set of test reads included in the inStrain source code.

### Step 1 — Download and prepare test data

Download the test files:

```bash
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1.R1.fastq.gz
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1.R2.fastq.gz
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1_scaffold_min1000.fa
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1.maxbin2.stb
curl -L -O https://raw.githubusercontent.com/MrOlm/inStrain/master/test/test_data/N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam
```

Or browse and download them manually from [GitHub](https://github.com/MrOlm/inStrain/tree/master/test/test_data).

| File | Description |
|------|-------------|
| `N5_271_010G1.R1.fastq.gz` | Forward metagenomic reads |
| `N5_271_010G1.R2.fastq.gz` | Reverse metagenomic reads |
| `N5_271_010G1_scaffold_min1000.fa` | Reference genome FASTA (from metagenomic assembly) |
| `N5_271_010G1.maxbin2.stb` | Scaffold-to-bin (STB) file mapping scaffolds to genomes |
| `N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam` | Pre-made BAM file for a second sample |

!!! info "About this sample"
    These reads come from a premature infant fecal sample.

Next, map the reads to the reference to generate a BAM file for the first sample:

```bash
mkdir bt2

bowtie2-build N5_271_010G1_scaffold_min1000.fa bt2/N5_271_010G1_scaffold_min1000.fa

bowtie2 -p 6 \
    -x bt2/N5_271_010G1_scaffold_min1000.fa \
    -1 N5_271_010G1.R1.fastq.gz \
    -2 N5_271_010G1.R2.fastq.gz \
    | samtools sort -@ 6 -o N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam
```

!!! warning "Dependencies"
    This step requires `bowtie2` and `samtools` to be installed and available on your `PATH`.

Finally, create a sample table listing both BAM files:

```bash
printf 'sample_name,bamfile\nN5_271_010G1,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam\nN5_271_010G2,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam\n' > samples.txt
```

The sample table should look like this:

```text
sample_name,bamfile
N5_271_010G1,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.bam
N5_271_010G2,N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.sorted.bam
```

---

### Step 2 — Profile

`zipstrain profile` needs only your BAM table, the reference FASTA, and the STB file. Any intermediate assets (null model, bed file, genome length table, profiling contract) are generated automatically into a `profiling_assets/` folder inside the run directory and reused on later runs:

```bash
zipstrain profile \
    --input-table samples.txt \
    --reference-fasta N5_271_010G1_scaffold_min1000.fa \
    --stb-file N5_271_010G1.maxbin2.stb \
    --run-dir out_profile
```

![zipstrain profile progress](assets/tutorial_profile.gif)

Results are written with one subdirectory per sample directly under the run directory, e.g. `out_profile/N5_271_010G1/`, each containing `*_profile.parquet`, `*_genome_stats.parquet`, and `*_gene_stats.parquet`. Per-sample intermediate files are tucked into an `intermediate_files/` subfolder, and the auto-generated profiling assets plus run logs live in `out_profile/profiling_assets/`.

!!! tip "Gene profiling and tuning"
    Pass `--gene-fasta` to enable gene-level profiling, or the null-model options (`--error-rate`, `--max-total-reads`, `--p-threshold`) to tune SNV calling. Use `--force-prepare` to rebuild the cached assets.

??? note "Preparing assets separately (advanced / cluster use)"
    You can still build the intermediate files ahead of time with `zipstrain utilities prepare_profiling` and pass them in explicitly via `--null-model`, `--bed-file`, `--genome-length-file`, etc. This is useful on clusters where you prepare once and profile many BAMs across nodes, and it is the path the Nextflow pipeline uses internally.

    ```bash
    zipstrain utilities prepare_profiling \
        --reference-fasta N5_271_010G1_scaffold_min1000.fa \
        -s N5_271_010G1.maxbin2.stb \
        -o profile_prep
    ```

---

### Step 3 — Compare

Build a profile database from the profile outputs:

```bash
printf 'profile_name,profile_location\nN5_271_010G1,out_profile/N5_271_010G1/N5_271_010G1_profile.parquet\nN5_271_010G2,out_profile/N5_271_010G2/N5_271_010G2_profile.parquet\n' > profiles.txt

zipstrain utilities build-profile-db \
    --profile-db-csv profiles.txt \
    --output-file profile_db.parquet
```

Then run the genome comparison:

```bash
zipstrain compare genomes \
    --profile-db profile_db.parquet \
    --scope all \
    --min-cov 5 \
    --stb-file N5_271_010G1.maxbin2.stb \
    --run-dir outputs/standard_compare \
    --calculate ani \
    --ani-method popani \
    --engine duckdb
```

Results are written to `outputs/standard_compare/Outputs/all_comparisons.parquet`. Each row is a genome-level pairwise comparison:

```text
Rows: 2
Columns: 6
$ genome           <str> 'maxbin2.maxbin.001.fasta', 'fobin.fasta'
$ total_positions  <i64> 0, 1183
$ share_allele_pos <i64> 0, 1183
$ genome_pop_ani   <f64> 0.0, 100.0
$ sample_1         <str> 'N5_271_010G1', 'N5_271_010G1'
$ sample_2         <str> 'N5_271_010G2', 'N5_271_010G2'
```

!!! info "Interpreting results"
    `genome_pop_ani` is the population-level ANI between the two samples for that genome. A value of `100.0` indicates the same strain is present in both samples; `0.0` indicates insufficient coverage to make a comparison.

TODO: Add an "expected results" section for interpretation of results

---

## Tutorial #2 — Nextflow with public SRA reads and auto-built genome database

This tutorial runs the same two samples as Tutorial #1, but pulls them directly from SRA and automatically builds a reference genome database using Sylph — no manual downloads or mapping required. This is the recommended route for larger studies or cluster environments.

!!! note "Prerequisites"
    Nextflow and Docker must be installed. See [Installation](./installation.md) for Nextflow setup. The `zipstrain.nf` pipeline file and `conf.config` are included in the [ZipStrain GitHub repo](https://github.com/OlmLab/ZipStrain).

### Step 1 — Create the SRA input table

The two samples used in this tutorial are available on NCBI SRA:

- [N5_271_010G1](https://www.ncbi.nlm.nih.gov/sra/?term=N5_271_010G1)
- [N5_271_010G2](https://www.ncbi.nlm.nih.gov/sra/?term=N5_271_010G2)

Find the **Run** accession (starts with `SRR`) for each on those pages, then create the input table:

```text
Run
SRR6262445
SRR6262448
```

```bash
printf 'Run\nSRR6262445\nSRR6262448\n' > sra.csv
```

---

### Step 2 — Profile (SRA to profile in one step)

This mode downloads reads from SRA, auto-builds a reference genome database via Sylph, maps reads, and generates profiles — all in one command:

```bash
nextflow run zipstrain.nf \
  --mode from_sra_to_profile \
  --input_table sra.csv \
  --output_dir out_sra_profile \
  -c conf.config \
  -profile docker \
  -resume
```

!!! tip "Auto-built genome database"
    Because no `--reference_genome` is provided, the pipeline runs Sylph on each sample to identify likely organisms and automatically downloads and builds the reference genome database. Results are cached in `--genome_db_cache_dir` (default: `genome_cache/`) so re-runs are fast.

Profiles are written to `out_sra_profile/profiles/`, one subdirectory per sample.

---

### Step 3 — Compare

Build a profile table pointing at the outputs from Step 2:

```bash
printf 'sample_name,profile_location\nN5_271_010G1,out_sra_profile/profiles/N5_271_010G1_profile.parquet\nN5_271_010G2,out_sra_profile/profiles/N5_271_010G2_profile.parquet\n' > profiles.csv
```

Then run genome comparison:

```bash
nextflow run zipstrain.nf \
  --mode compare_genomes \
  --input_type profile_table \
  --input_table profiles.csv \
  --compare_genome_scope all \
  --compare_calculate ani \
  --compare_ani_method popani \
  --compare_engine duckdb \
  --output_dir out_sra_compare \
  -c conf.config \
  -profile docker \
  -resume
```

Results are written to `out_sra_compare/Outputs/all_comparisons.parquet` in the same format as Tutorial #1.
