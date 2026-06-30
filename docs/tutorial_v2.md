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

First, prepare the assets needed for profiling:

```bash
zipstrain utilities prepare_profiling \
    --reference-fasta N5_271_010G1_scaffold_min1000.fa \
    -s N5_271_010G1.maxbin2.stb \
    -o profile_prep
```

??? tip "Additional options for `prepare_profiling`"
    ```bash
    zipstrain utilities prepare_profiling --help
    ```

    ```text
    Usage: zipstrain utilities prepare_profiling [OPTIONS]

      Prepare the files needed for profiling bam files and save them in the
      specified output directory.

    Options:
      -r, --reference-fasta TEXT     Path to the reference genome in FASTA format.
                                     [required]
      -g, --gene-fasta TEXT          Optional path to the gene annotations in
                                     FASTA format.
      -s, --stb-file TEXT            Path to the scaffold-to-genome mapping file.
                                     [required]
      -e, --error-rate FLOAT         Error rate for the sequencing technology when
                                     building the null model.  [default: 0.001]
      -m, --max-total-reads INTEGER  Maximum coverage to consider when building
                                     the null model.  [default: 10000]
      -p, --p-threshold FLOAT        Significance threshold for the null model.
                                     [default: 0.05]
      -t, --model-type [poisson]     Type of null model to build.  [default:
                                     poisson]
      -o, --output-dir TEXT          Directory to save the profiling database.
                                     [required]
      --help                         Show this message and exit.
    ```

!!! note "Gene profiling"
    If you'd like to profile genes or adjust SNV-calling parameters, pass `--gene-fasta` to `prepare_profiling` before running `profile`.

Then run `zipstrain profile`:

```bash
zipstrain profile \
    --input-table samples.txt \
    --reference-fasta N5_271_010G1_scaffold_min1000.fa \
    --stb-file N5_271_010G1.maxbin2.stb \
    --null-model profile_prep/null_model.parquet \
    --profiling-contract profile_prep/profiling_contract.json \
    --bed-file profile_prep/genomes_bed_file.bed \
    --genome-length-file profile_prep/genome_lengths.parquet \
    --run-dir out_profile
```

![zipstrain profile progress](assets/tutorial_profile.gif)

Results will be written to `out_profile/batch_0/`, with one subdirectory per sample containing a `*_genome_stats.parquet` and `*_profile.parquet`.

---

### Step 3 — Compare

Build a profile database from the profile outputs:

```bash
printf 'profile_name,profile_location\nN5_271_010G1,out_profile/batch_0/N5_271_010G1/N5_271_010G1_profile.parquet\nN5_271_010G2,out_profile/batch_0/N5_271_010G2/N5_271_010G2_profile.parquet\n' > profiles.txt

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
<SRR_accession_for_N5_271_010G1>
<SRR_accession_for_N5_271_010G2>
```

```bash
printf 'Run\n<SRR_G1>\n<SRR_G2>\n' > sra.csv
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
