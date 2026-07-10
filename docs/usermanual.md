# User Manual

##  An overview of ZipStrain and the data it generates

ZipStrain is a program for microbial metagenomic analysis. The primary use cases for ZipStrain are to accurately determine organism presence / absence in a community, and to perform detailed comparisons between organisms in different samples.

There are two ways of interacting with ZipStrain: the **Python CLI** and the **Nextflow pipeline**. Both run the same underlying analysis and produce the same output tables; they differ in how the work is orchestrated and how much infrastructure they assume.

**Python CLI** (`zipstrain map | profile | compare`)

- *Pros:* nothing to install beyond ZipStrain and its dependencies; easy to run one step at a time and inspect intermediates; simplest path on a laptop, a single workstation, or an interactive session; the newer matrix-store comparison workflow is CLI-only.
- *Cons:* you drive the steps yourself; large cohorts on a cluster mean managing SLURM submission and resumption by hand (the commands support `--execution-mode slurm`, but Nextflow does more of this for you).
- *Best when:* you are exploring, working on a handful of samples, prototyping, or want the matrix workflow.

**Nextflow pipeline** (`nextflow run zipstrain.nf`)

- *Pros:* one command runs the whole map → profile → compare chain; built-in `-resume`, containerization (Docker/Singularity/Apptainer), and scheduler execution (SLURM, etc.) via execution profiles; scales cleanly to large cohorts and shared HPC systems; reproducible across machines.
- *Cons:* requires Nextflow (and usually a container engine) and some config setup; less convenient for poking at a single intermediate step; the matrix comparison workflow is not yet wired into Nextflow.
- *Best when:* you are processing many samples, running on a cluster, or want a reproducible, restartable, containerized pipeline.

A common pattern is to prototype with the CLI on a few samples, then scale the same analysis out with Nextflow. The [ZipStrain Command Line Interface](#zipstrain-command-line-interface) section below documents the CLI; [Nextflow Implementation](#nextflow-implementation) covers the pipeline.

TODO: Make a figure like this inStrain one (or just edit this inStrain one) so people have a vauge idea of what's going on from the bat: https://instrain.readthedocs.io/en/latest/_images/OverviewFigure1_v1.4.png . We can then add this figure lots of places

## ZipStrain Command Line Interface

A typical ZipStrain run has three steps, each its own top-level command:

1. **`zipstrain map`** — turn reads into sorted BAM files
2. **`zipstrain profile`** — profile those BAMs at nucleotide resolution
3. **`zipstrain compare`** — compare samples to each other by ANI

Two more commands round it out: **`zipstrain test`** checks your environment, and **`zipstrain utilities`** groups the lower-level helpers that the three main commands are built from (you rarely call these directly).

For end-to-end, copy-pasteable examples, start with the [Tutorial](./tutorial_v2.md); for the files each command writes, see [Expected output](./expected_output.md). This page is the reference for what each command does and its major options.

Every command supports `-h`/`--help`, and grouped utilities take a subcommand:

```bash
zipstrain --help                 # top-level commands
zipstrain <command> --help       # e.g. zipstrain profile --help
zipstrain utilities <cmd> --help # e.g. zipstrain utilities build-genome-db --help
```

<details>
<summary><code>zipstrain --help</code></summary>

```text
Usage: zipstrain [OPTIONS] COMMAND [ARGS]...

  ZipStrain — fast strain-level metagenomic profiling and comparison.

  A typical run goes: map reads to BAMs, profile them at nucleotide
  resolution, then compare samples by ANI.

  Developed by Parsa Ghadermazi and Matt Olm in the Olm Lab at the University
  of Colorado Boulder.

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  test       Check your environment is ready.
  map        Map reads to sorted BAMs.
  profile    Profile BAMs at nucleotide resolution.
  compare    Compare samples (genomes or genes).
  utilities  Lower-level helper commands.

  Source & docs: https://github.com/OlmLab/ZipStrain
```

</details>

### Run logging

Every `map`, `profile`, and `compare` run writes two logging files into its output directory (`--output-dir` for `map`, `--run-dir` for `profile`/`compare`), so you can tell at a glance whether a run is still going, finished, or crashed:

- **`zipstrain_run.log`** — a human-readable, append-only timeline: a header (ZipStrain version, full command, working directory, process id, start time) followed by a timestamped line for each step. On success it ends with `COMPLETED`; on failure it ends with `CRASHED` (or `ABORTED` for a Ctrl-C), the step that was running, and the full Python traceback.
- **`zipstrain_status.json`** — a machine-readable snapshot with `status` (`running` / `completed` / `crashed` / `aborted`), `last_step`, `pid`, timestamps, `elapsed_seconds`, and — on failure — `error_type`, `error_message`, and `failed_during_step`.

A quick way to check on a run:

```bash
# What is (or was) it doing?
grep -E 'RUNNING|STEP|COMPLETED|CRASHED|ABORTED' out_dir/zipstrain_run.log | tail
# Just the status
python -c "import json; print(json.load(open('out_dir/zipstrain_status.json'))['status'])"
```

If `status` is `running` but nothing has updated for a long time, use the recorded `pid` to check whether the process is still alive (e.g. `ps -p <pid>`). Because the log is append-only, re-running in the same directory (for example to [resume `map`](#map)) adds a fresh block rather than overwriting the previous run's history.

### Map

#### `zipstrain map`

Turn sequencing reads into sorted, indexed BAM files ready for `zipstrain profile`. You provide a reads table (`sample_name,reads1[,reads2]`). If you also pass `--reference-fasta` (plus its `--stb-file`), ZipStrain maps against it with Bowtie2. If you **omit** the reference, ZipStrain runs Sylph to auto-pick reference genomes from the reads, downloads and caches them, and maps against the built reference — so a minimal run needs only a reads table, an output directory, and (for the Sylph route) a database and cache directory.

Alongside the BAMs it writes the reference FASTA + STB and a `samples.txt` you can hand straight to `zipstrain profile`. On the Sylph route it also writes a `reference_genomes_taxonomy.tsv` that `profile` auto-discovers to populate `genome_taxonomy`.

```bash
# Sylph auto-reference (no reference on hand)
zipstrain map \
  --reads-table reads.csv \
  --output-dir out_map \
  --sylph-db gtdb-r220-c200-dbv1.syldb \
  --genome-cache-dir genome_cache

# Map against a reference you already have
zipstrain map \
  --reads-table reads.csv \
  --reference-fasta reference_genomes.fna \
  --stb-file reference_genomes.stb \
  --output-dir out_map
```

Major options:

- `-i, --reads-table` / `-o, --output-dir` (required) — reads CSV and output directory
- `-f, --reference-fasta` + `-s, --stb-file` — map against an existing reference (skips Sylph)
- `--sylph-db` / `--sylph-db-url` / `--genome-cache-dir` — Sylph database and genome cache for the auto-reference route
- `--predict-genes` — also run prodigal to emit a gene FASTA for later gene-level profiling
- `--non-competitive` — pass `-a` to Bowtie2 (report all alignments)
- `--force` — redo every step even if cached outputs exist (see resume note below)
- `-t, --threads` (default: `4`)

`zipstrain map` is **resumable**: re-running with the same `--output-dir` skips any stage whose output is already complete — per-sample Sylph tables, the built reference, the Bowtie2 index, and each sorted+indexed BAM. So if a run crashes partway through (say, on the third sample), just run the same command again and it picks up where it left off. Completion is judged conservatively (a BAM only counts as done once its `.bai` index exists, writes are atomic), so a half-written file from a crash is never reused. Pass `--force` to ignore the cache and rebuild everything.

<details>
<summary><code>zipstrain map --help</code></summary>

```text
Usage: zipstrain map [OPTIONS]

  Map sequencing reads to BAM files, ready for `zipstrain profile`.

  Provide a reads table (``sample_name,reads1[,reads2]``). If you do not pass
  ``--reference-fasta``, ZipStrain runs Sylph to pick reference genomes from
  the reads automatically, downloading and caching them, then maps against the
  built reference. Outputs sorted, indexed BAMs, the reference FASTA + STB,
  and a ``samples.txt`` you can hand straight to ``zipstrain profile``.

Required inputs:
  -i, --reads-table TEXT  CSV of reads to map, with columns
                          'sample_name,reads1[,reads2]' (reads2 blank/absent
                          for single-end).  [required]
  -o, --output-dir TEXT   Directory to write BAMs, the reference FASTA/STB,
                          and a samples.txt ready for `zipstrain profile`.
                          [required]

Reference (Sylph auto-picks one if omitted):
  -f, --reference-fasta TEXT  Reference FASTA to map against. If omitted,
                              Sylph automatically picks and builds a reference
                              from the reads.
  -s, --stb-file TEXT         Scaffold-to-genome mapping file. Required when
                              --reference-fasta is provided.
  --sylph-db TEXT             Path to the Sylph database. Used when no
                              --reference-fasta is given; downloaded from
                              --sylph-db-url if the path does not exist.
  --sylph-db-url TEXT         URL to download the Sylph database from when
                              --sylph-db is missing.  [default:
                              http://faust.compbio.cs.cmu.edu/sylph-
                              stuff/gtdb-r220-c200-dbv1.syldb]
  --genome-cache-dir TEXT     Directory that caches genome FASTAs downloaded
                              during Sylph-based reference building. Required
                              when no --reference-fasta is given.

Options:
  --predict-genes        Also run prodigal to emit a gene FASTA (for gene-
                         level profiling via `profile --gene-fasta`).
  --non-competitive      Pass -a to Bowtie2 for non-competitive mapping
                         (report all alignments).
  --force                Redo every step from scratch, ignoring cached
                         outputs. By default `map` resumes: completed Sylph
                         tables, reference, index, and BAMs are reused.
  -t, --threads INTEGER  Threads for Sylph, Bowtie2, and samtools.  [default:
                         4]

Other options:
  -h, --help  Show this message and exit.
```

</details>

### Profile

#### `zipstrain profile`

Profile a batch of BAM files at nucleotide resolution, producing per-position base counts plus per-genome and per-gene summary tables. Profiling needs a set of assets (null model, BED file, genome-length table, optional gene ranges, profiling contract). You can supply any of these, but any you omit are **generated automatically** into a `profiling_assets/` directory inside `--run-dir` and reused on later runs when the inputs are unchanged. A minimal run therefore needs only `--input-table`, `--reference-fasta`, and `--stb-file`.

By default each sample also gets an SNV table (`<sample>_SNVs.parquet`) and a `presence` (present/absent) call in its genome stats. Companion `.csv` files are written next to the small stat tables. See [Expected output](./expected_output.md) for the full file list.

```bash
# Minimal run — assets auto-generated and cached in run_dir/profiling_assets
zipstrain profile \
  --input-table samples.txt \
  --reference-fasta reference_genomes.fna \
  --stb-file reference_genomes.stb \
  --run-dir out_profile
```

Major options:

- `-i, --input-table`, `-f, --reference-fasta`, `-s, --stb-file`, `-r, --run-dir` (required)
- `--gene-fasta` — enables gene-level profiling (auto-generates a gene range table)
- `-u/-b/-l/-g/--profiling-contract` — supply pre-built assets to override auto-generation; `--force-prepare` regenerates them all
- Read filters: `--min-mapq` (0), `--min-baseq` (13), `--min-read-ani`, `--read-inclusion` (`all-mapped`)
- Presence & SNVs: `--no-snvs`, `--snv-min-cov` (5), `--presence-ber` (0.5), `--presence-fug` (1.0), `--presence-min-cov-use-fug` (2.0), `--presence-min-coverage` (0.1), `--genome-taxonomy`
- Output: `--no-csv` / `--force-csv` (companion CSVs; 100 MB cap by default)
- Execution: `-n, --num-procs` (8), `-m, --max-concurrent-batches` (5), `-t, --task-per-batch` (10), `-e, --execution-mode` (`local`/`slurm`), `-c, --slurm-config`, `-o, --container-engine`, `--container-address`

<details>
<summary><code>zipstrain profile --help</code></summary>

```text
Usage: zipstrain profile [OPTIONS]

  Run BAM file profiling in batches using the specified execution mode and
  container engine.

  Any profiling assets (null model, bed file, genome length table, gene range
  table, profiling contract) that are not supplied explicitly are generated
  automatically into a ``profiling_assets`` directory inside ``run-dir`` and
  reused on subsequent runs when the inputs are unchanged. This means a
  minimal run needs only ``--input-table``, ``--reference-fasta``, and
  ``--stb-file``.

Required inputs:
  -i, --input-table TEXT      Path to the input table in TSV format containing
                              sample names and paths to bam files.  [required]
  -f, --reference-fasta TEXT  Reference FASTA. Used for mpileup (adds
                              ref_base_bitmask) and required to auto-generate
                              the bed/genome-length assets when they are not
                              supplied.
  -s, --stb-file TEXT         Path to the scaffold-to-genome mapping file.
                              [required]
  -r, --run-dir TEXT          Directory to save the run data (sample outputs,
                              profiling_assets, and logs).  [required]

Optional inputs:
  --gene-fasta TEXT  Gene FASTA. When provided, a gene range table is auto-
                     generated from it for gene-level profiling.

Optional pre-built assets (auto-generated if omitted):
  -u, --null-model TEXT          Pre-built null model parquet file. Auto-
                                 generated into <run-dir>/profiling_assets if
                                 not provided.
  -g, --gene-range-table TEXT    Pre-built gene range table file. Overrides
                                 --gene-fasta auto-generation.
  --profiling-contract TEXT      Pre-built profiling_contract.json. When
                                 provided, its hashes are written into each
                                 profile parquet metadata. Auto-generated
                                 otherwise.
  -b, --bed-file TEXT            Pre-built BED file for profiling regions.
                                 Auto-generated into <run-
                                 dir>/profiling_assets if not provided.
  -l, --genome-length-file TEXT  Pre-built genome length file. Auto-generated
                                 into <run-dir>/profiling_assets if not
                                 provided.

Profiling parameters:
  --error-rate FLOAT              Error rate used when auto-generating the
                                  null model.  [default: 0.001]
  --max-total-reads INTEGER       Maximum coverage considered when auto-
                                  generating the null model.  [default: 10000]
  --p-threshold FLOAT             Significance threshold used when auto-
                                  generating the null model.  [default: 0.05]
  --model-type [poisson]          Null model type used when auto-generating
                                  the null model.  [default: poisson]
  --force-prepare                 Regenerate all auto-generated profiling
                                  assets even if valid cached copies exist.
  --min-mapq INTEGER              Minimum mapping quality for a read to be
                                  used during profiling.  [default: 0]
  --min-baseq INTEGER             Minimum base quality for a base to be
                                  counted during profiling.  [default: 13]
  --min-read-ani FLOAT            Minimum read ANI proxy based on the NM tag
                                  and aligned query span.
  --read-inclusion [proper-pairs|paired|all-mapped]
                                  Which mapped reads are eligible for
                                  profiling.  [default: all-mapped]

SNV calling and presence:
  --no-snvs                       Do not call SNVs/SNPs (per-sample
                                  <sample>_SNVs.parquet). SNV calling needs
                                  --reference-fasta.
  --snv-min-cov INTEGER           Minimum coverage for a site to be eligible
                                  as an SNV/SNP call.  [default: 5]
  --presence-ber FLOAT            Breadth-error-ratio threshold for the genome
                                  present/absent call (the Metapresence paper
                                  recommends ~0.8).  [default: 0.5]
  --presence-fug FLOAT            FUG threshold for the present/absent call at
                                  low coverage. A genome is present when
                                  fug/0.632 exceeds this (fug ~ 0.632 under
                                  uniform coverage, so 1.0 means at least as
                                  uniform as random).  [default: 1.0]
  --presence-min-cov-use-fug FLOAT
                                  Coverage above which the present/absent call
                                  uses BER alone (below it, FUG is also
                                  required).  [default: 2.0]
  --presence-min-coverage FLOAT   Minimum mean coverage required to call a
                                  genome present.  [default: 0.1]
  --genome-taxonomy TEXT          Optional genome->taxonomy TSV to add a
                                  genome_taxonomy column to genome_stats.
                                  Auto-discovered next to the reference/STB
                                  when produced by `zipstrain map` (Sylph
                                  route).

Output:
  --no-csv     Do not write companion .csv files next to the
               genome_stats/gene_stats/SNV parquets.
  --force-csv  Write companion .csv files even when the estimated size exceeds
               100 MB.

Running parameters:
  -n, --num-procs INTEGER         Number of processors to use for each
                                  profiling task.  [default: 8]
  -m, --max-concurrent-batches INTEGER
                                  Maximum number of concurrent batches to run.
                                  [default: 5]
  -p, --poll-interval INTEGER     Polling interval in seconds to check the
                                  status of batches.  [default: 1]
  -t, --task-per-batch INTEGER    Number of tasks to include in each batch.
                                  [default: 10]
  -e, --execution-mode TEXT       Execution mode: 'local' or 'slurm'.
                                  [default: local]
  -c, --slurm-config TEXT         Path to the SLURM configuration file in json
                                  format. Required if execution mode is
                                  'slurm'.
  -o, --container-engine TEXT     Container engine to use: 'local', 'docker'
                                  or 'apptainer'.  [default: local]
  --container-address TEXT        Optional container image/address override.
                                  Defaults to the current ZipStrain version
                                  tag for docker/apptainer.

Other options:
  -h, --help  Show this message and exit.
```

</details>

### Comparison

#### `zipstrain compare`

Compare profiled samples to each other, one row per genome per sample pair, and write `<run-dir>/all_comparisons.parquet` (+ a companion CSV). By default it compares at the genome level; add `--compare-genes` for gene-level comparison. `--profile-db` accepts a CSV of `profile_name,profile_location` rows directly (no need to run `build-profile-db` first) or a pre-built profile-database parquet.

There are two engines. `--method standard` (default) does direct pairwise comparison and is simplest. `--method matrix` builds a reusable matrix store, which pays off for repeated all-vs-all comparison. Both are **resumable and extendable**: re-running with the same `--run-dir` and a profiles table that includes new samples computes only the new pairs. See the [Tutorial](./tutorial_v2.md) for a worked matrix walk-through and [Expected output](./expected_output.md) for the columns.

```bash
# Standard genome comparison from a CSV of profiles
zipstrain compare \
  --profile-db profiles.csv \
  --run-dir out_compare

# Matrix method, reusable for repeated all-vs-all
zipstrain compare \
  --profile-db profiles.csv \
  --run-dir out_compare \
  --method matrix \
  --stb-file reference_genomes.stb
```

Major options:

- `--profile-db`, `-r, --run-dir` (required)
- `--method` (`standard`/`matrix`) and `--compare-genes` — pick the engine and genome vs. gene level
- `--scope` (`all` for genomes, `all:all` for genes), `--min-cov` (5), `--min-gene-compare-len` (100)
- `-a, --ani-method` (`popani`, `conani`, `cosani_<threshold>`), `--calculate` (`ani`/`ibs`/`identical_genes`/`all`)
- `--stb-file` — required for `--method matrix`
- Matrix method: `--bed-file`, `--gene-range-table` (both auto-discovered from `profiling_assets`), `--backend` (`numpy`/`torch…`), `--memory-limit-gb` (16)
- Standard method: `--engine` (`polars`/`duckdb`), `-d, --duckdb-memory-limit`, `--duckdb-threads`, plus the same execution/container options as `profile`
- Output: `--no-csv` / `--force-csv`
- `--comp-db-file` / `--allow-mismatch` — resume from an existing comparison / skip profile-contract validation

<details>
<summary><code>zipstrain compare --help</code></summary>

```text
Usage: zipstrain compare [OPTIONS]

  Compare profiled samples at the genome level (default) or gene level
  (--compare-genes).

  ``--profile-db`` may be a CSV of ``profile_name,profile_location`` rows, so
  there is no need to run ``zipstrain utilities build-profile-db`` first; a
  pre-built profile-database parquet is also accepted.

  Both methods write ``<run-dir>/all_comparisons.parquet``. Re-running with
  the same ``--run-dir`` and a profiles table that includes new samples
  extends the existing comparison, computing only the new pairs.

Required inputs:
  --profile-db TEXT   Profiles to compare: either a CSV with
                      'profile_name,profile_location' columns (built in
                      memory, no build-profile-db needed) or a pre-built
                      profile-database parquet.  [required]
  -r, --run-dir TEXT  Directory to save the run data.  [required]

Comparison parameters:
  --method [standard|matrix]      Comparison engine: 'standard' (direct
                                  pairwise) or 'matrix' (reusable matrix
                                  store, good for repeated all-vs-all).
                                  [default: standard]
  --compare-genes                 Compare genes instead of genomes.
  --scope TEXT                    Comparison scope. Defaults to 'all' for
                                  genomes and 'all:all' for genes.
  --min-cov INTEGER               Minimum coverage to consider a position.
                                  [default: 5]
  --min-gene-compare-len INTEGER  Minimum gene length to consider for
                                  comparison.  [default: 100]
  --stb-file TEXT                 Scaffold-to-genome mapping file. Required
                                  for --method matrix.
  -a, --ani-method TEXT           ANI calculation method (e.g., 'popani',
                                  'conani', 'cosani_0.4').  [default: popani]
  --calculate TEXT                Genome metrics to compute (genome mode
                                  only): ani, ibs, identical_genes. Combine
                                  with '+', or use all.  [default: all]
  --comp-db-file TEXT             Optional existing comparison parquet to
                                  resume/extend (standard method). Auto-
                                  detected from --run-dir if omitted.
  --allow-mismatch                Skip profile contract validation when
                                  building the profile database from a CSV.

Matrix method (--method matrix):
  --bed-file TEXT                 BED file for the matrix store (--method
                                  matrix). Auto-discovered from
                                  profiling_assets if omitted.
  --gene-range-table TEXT         Gene range table for gene ANI (--method
                                  matrix). Auto-discovered from
                                  profiling_assets if omitted.
  --backend [numpy|torch|torch-cpu|torch-cuda|torch-mps]
                                  Compute backend for --method matrix (numpy,
                                  or torch on CPU/CUDA/MPS).  [default: numpy]
  --memory-limit-gb FLOAT         Approximate memory budget for --method
                                  matrix.  [default: 16.0]

Output:
  --no-csv     Do not write a companion .csv next to the comparison parquet.
  --force-csv  Write the companion .csv even when the estimated size exceeds
               100 MB.

Standard method / engine:
  --engine [polars|duckdb]        Comparison engine for standard compare
                                  tasks.  [default: polars]
  -d, --duckdb-memory-limit TEXT  DuckDB memory limit for compare tasks (e.g.,
                                  2GB).
  --duckdb-threads INTEGER        Number of DuckDB worker threads for compare
                                  tasks.
  -m, --max-concurrent-batches INTEGER
                                  Maximum number of concurrent batches to run.
                                  [default: 5]
  -p, --poll-interval INTEGER     Polling interval in seconds to check the
                                  status of batches.  [default: 1]
  -t, --task-per-batch INTEGER    Number of tasks to include in each batch.
                                  [default: 10]
  -e, --execution-mode TEXT       Execution mode: 'local' or 'slurm'.
                                  [default: local]
  -s, --slurm-config TEXT         Path to the SLURM configuration file in json
                                  format. Required if execution mode is
                                  'slurm'.
  -c, --container-engine TEXT     Container engine to use: 'local', 'docker'
                                  or 'apptainer'.  [default: local]
  --container-address TEXT        Optional container image/address override.
                                  Defaults to the current ZipStrain version
                                  tag for docker/apptainer.

Other options:
  -h, --help  Show this message and exit.
```

</details>

### Test

#### `zipstrain test`

Run a lightweight health check that confirms ZipStrain and its external dependencies are importable and callable. Run it right after installation. It takes no options.

```bash
zipstrain test
```

<details>
<summary><code>zipstrain test --help</code></summary>

```text
Usage: zipstrain test [OPTIONS]

  Run a lightweight ZipStrain health check.

Options:
  -h, --help  Show this message and exit.
```

</details>

### Utilities

The `zipstrain utilities` group collects the lower-level helpers that `map`, `profile`, and `compare` are built from — asset preparation, single-BAM/single-pair operations, matrix-store management, format conversions, and more. Most users never call these directly, but they are useful for building custom pipelines or running one step in isolation. Use `zipstrain utilities <command> --help` for full details of any one.

<details>
<summary><code>zipstrain utilities --help</code></summary>

```text
Usage: zipstrain utilities [OPTIONS] COMMAND [ARGS]...

  The commands in this group are related to various utility functions that
  mainly prepare input files for profiling and comparison.

Options:
  -h, --help  Show this message and exit.

Commands:
  adjust-sequence-errors  Apply ZipStrain's sequence-error adjustment to...
  append-matrix-db        Append new profiles to an existing matrix store.
  build-genome-db         Build a reference bundle from an abundance table.
  build-matrix-db         Build a matrix store directly from classic...
  build-null-model        Build a null model for sequencing errors based...
  build-profile-db        Build a profile database from the given CSV file.
  chunk-genome-compare    Run classic genome compare over one pair-table...
  gene-range-table        Main function to build and save the gene...
  generate-genome-pairs   Generate a pair table ready for...
  generate_stb            Generate a scaffold-to-genome mapping file from...
  get-coverage-stats      Build coverage-only gene and genome stats from...
  get-snp-reference       Emit profile-like rows that are SNPs relative...
  get_genome_lengths      Extract the genome length information from the...
  make_bed                Create a BED file from the database in fasta...
  matrix-compare          Run resumable matrix compare on all...
  matrix-compare-export   Export a matrix compare DuckDB database to...
  matrix-db-to-hdf5       Convert a legacy DuckDB matrix database into...
  merge-stat-tables       Concatenate stat tables and add a sample column...
  merge_parquet           Merge multiple Parquet files in a directory...
  prepare_profiling       Prepare the files needed for profiling bam...
  presence-profile        Generate a presence profile for genomes based...
  process-read-locs       Process read locations and save them to a...
  process_mpileup         Process mpileup files and save the results in a...
  profile-single          Profile a single BAM file using the provided...
  single_compare_gene     Compare two mpileup files and calculate...
  single_compare_genome   Main function to compare two mpileup files and...
  sort-profile            Sort a classic profile parquet in place and...
  strain_heterogeneity    Calculate strain heterogeneity for each genome...
  to-complete-table       Generate the not-yet-completed...
```

</details>

**Utility Commands At A Glance**

| Command | Purpose |
|---|---|
| `zipstrain utilities build-null-model` | Build sequencing-error null model |
| `zipstrain utilities merge_parquet` | Merge parquet files |
| `zipstrain utilities merge-stat-tables` | Merge gene/genome stat parquet files with sample labels |
| `zipstrain utilities get-coverage-stats` | Rebuild coverage-only gene/genome stats from a profile parquet |
| `zipstrain utilities process_mpileup` | Convert mpileup stream to parquet |
| `zipstrain utilities make_bed` | Build bed chunks from fasta |
| `zipstrain utilities get_genome_lengths` | Genome lengths from STB + BED |
| `zipstrain utilities generate-genome-pairs` | Create all non-redundant standard-profile pairs |
| `zipstrain utilities chunk-genome-compare` | Compare many genome-level pairs in Python-side parallel batches |
| `zipstrain utilities strain_heterogeneity` | Strain heterogeneity metrics |
| `zipstrain utilities build-profile-db` | Build profile DB parquet |
| `zipstrain utilities build-matrix-db` | Build the current per-sample genome matrix store directly from profile parquets |
| `zipstrain utilities append-matrix-db` | Append new profiles into an existing matrix store |
| `zipstrain utilities matrix-db-to-hdf5` | Convert a DuckDB matrix database into the current matrix-store format |
| `zipstrain utilities matrix-compare` | Resumable all-vs-all matrix compare into a DuckDB compare DB |
| `zipstrain utilities matrix-compare-export` | Export a matrix compare DuckDB to parquet |
| `zipstrain utilities build-genome-db` | Build local genome reference bundle from abundance table |
| `zipstrain utilities presence-profile` | Presence profile from coverage + read locations |
| `zipstrain utilities process-read-locs` | Process read-location stream |
| `zipstrain utilities generate_stb` | Create scaffold-to-genome map from genome files |
| `zipstrain utilities gene-range-table` | Create gene range table |
| `zipstrain test` | Validate local installation/dependencies |

The commands below run individual profiling and comparison steps by hand (the `map`, `profile`, and `compare` commands orchestrate them for you).

#### `zipstrain utilities prepare_profiling`

Prepare profiling database assets.

```bash
zipstrain utilities prepare_profiling \
  --reference-fasta reference.fasta \
  --stb-file mapping.stb \
  --output-dir profiling_assets
```

Options:

- `-r, --reference-fasta` (required)
- `-g, --gene-fasta` (optional)
- `-s, --stb-file` (required)
- `-e, --error-rate` (default: `0.001`)
- `-m, --max-total-reads` (default: `10000`)
- `-p, --p-threshold` (default: `0.05`)
- `-t, --model-type` (default: `poisson`)
- `-o, --output-dir` (required)

Outputs:

- `genomes_bed_file.bed`
- `gene_range_table.tsv`
- `genome_lengths.parquet`
- `null_model.parquet`
- `profiling_contract.json`

#### `zipstrain utilities profile-single`

Profile a single BAM.

```bash
zipstrain utilities profile-single \
  --bed-file genomes_bed_file.bed \
  --bam-file sample.bam \
  --stb-file mapping.stb \
  --null-model profiling_assets/null_model.parquet \
  --profiling-contract profiling_assets/profiling_contract.json \
  --num-chunks 24 \
  --max-concurrency 4 \
  --output-dir sample_profile
```

Options:

- `-r, --reference-fasta` (optional) — when provided, profiling also records `ref_base_bitmask` and adds `ref_ani` to gene/genome stat tables
- `-b, --bed-file` (required)
- `-a, --bam-file` (required)
- `-s, --stb-file` (required)
- `-m, --null-model` (required)
- `-g, --gene-range-table` (optional)
- `--profiling-contract` (optional)
- `-n, --num-chunks` (default: `24`) — number of BED chunks to create
- `-c, --max-concurrency` (default: `4`) — how many chunks run simultaneously
- `--min-mapq` (default: `0`)
- `--min-baseq` (default: `13`)
- `--min-read-ani` (optional) — filters reads before pileup using the BAM `NM` tag and aligned query span
- `--read-inclusion` (`proper-pairs|paired|all-mapped`, default: `all-mapped`)
- `-o, --output-dir` (required)

Read-inclusion modes:

- `proper-pairs`: keep only mapped read pairs carrying the aligner `PROPER_PAIR` flag
- `paired`: keep mapped paired-end reads even if they are discordant
- `all-mapped`: keep any mapped read, whether paired or single-end

Outputs include:

- `<sample>_profile.parquet`
- `<sample>_genome_stats.parquet`
- `<sample>_gene_stats.parquet`

When `--reference-fasta` is provided during profiling, the profile parquet includes `ref_base_bitmask`.
In the same case, the generated genome and gene stat tables also include a `ref_ani` column.

`ref_ani` is the percentage of covered sites whose observed allele set still contains the reference allele after ZipStrain's sequence-error adjustment.

`ref_base_bitmask` uses this encoding:

- `1` = reference base `A`
- `2` = reference base `C`
- `4` = reference base `G`
- `8` = reference base `T`
- `0` = non-ACGT or unknown reference base

This is a one-hot bitmask, so current profiles are expected to contain only `0`, `1`, `2`, `4`, or `8` in this column.

#### `zipstrain utilities get-snp-reference`

Emit profile-like rows that are SNPs relative to the reference from a classic profile parquet that includes `ref_base_bitmask`.

```bash
zipstrain utilities get-snp-reference \
  --profile-file sample_profile.parquet \
  --min-cov 5 \
  --output-file sample_reference_snps.parquet
```

Options:

- `-p, --profile-file` (required)
- `-c, --min-cov` (default: `5`)
- `-o, --output-file` (required)

The output preserves the input profile-like columns and includes only positions that:

- have coverage `>= min_cov`
- have a known reference base in `ref_base_bitmask`
- do not retain the reference allele after profile sequence-error adjustment

This uses the same reference-sharing logic used to populate `ref_ani` in the gene and genome stat tables.

#### `zipstrain utilities single_compare_genome`

```bash
zipstrain utilities single_compare_genome \
  --mpileup-contig-1 sample_a.parquet \
  --mpileup-contig-2 sample_b.parquet \
  --stb-file mapping.stb \
  --output-file out.parquet
```

Options:

- `-m1, --mpileup-contig-1` (required)
- `-m2, --mpileup-contig-2` (required)
- `-s, --stb-file` (required)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `100`)
- `-o, --output-file` (required)
- `-g, --genome` (default: `all`)
- `-a, --ani-method` (default: `popani`)
- `--calculate` (default: `all`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`

#### `zipstrain utilities generate-genome-pairs`

```bash
zipstrain utilities generate-genome-pairs \
  --profile-dir profiles \
  --output-file genome_pairs.parquet
```

This writes a parquet table with:

- `sample_name_1`
- `sample_name_2`
- `profile_location_1`
- `profile_location_2`

Options:

- `-p, --profile-dir` (required)
- `-o, --output-file` (required)
- `--write-batch-size` (default: `100000`)

#### `zipstrain utilities chunk-genome-compare`

```bash
zipstrain utilities chunk-genome-compare \
  --pair-table genome_pairs.parquet \
  --stb-file mapping.stb \
  --output-file chunk_compare.parquet \
  --workers 8 \
  --engine polars
```

This command runs standard genome comparisons directly inside Python for one pair-table chunk.
It is intended as an experimental utility for benchmarking or ad hoc compare runs, and
does not change the main workflow commands.

Accepted pair-table schemas:

- `sample_name_1`, `sample_name_2`, `profile_location_1`, `profile_location_2`
- `sample_name_1`, `sample_name_2`, `profile_1`, `profile_2`
- `sample_1`, `sample_2`, `profile_1`, `profile_2`
- `profile_location_1`, `profile_location_2`
- `profile_1`, `profile_2`

Options:

- `-p, --pair-table` (required)
- `-s, --stb-file` (required)
- `-o, --output-file` (required)
- `-w, --workers` (defaults to CPU count capped by pair count)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `100`)
- `-g, --genome` (default: `all`)
- `-a, --ani-method` (default: `popani`)
- `--calculate` (default: `all`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`

The final console summary includes:

- total pairs processed
- total genome-level output rows written
- total elapsed time
- average wall time per pair
- average compute time per pair
- average time per genome-level output row

#### `zipstrain utilities single_compare_gene`

```bash
zipstrain utilities single_compare_gene \
  --mpileup-contig-1 sample_a.parquet \
  --mpileup-contig-2 sample_b.parquet \
  --stb-file mapping.stb \
  --scope all:all \
  --output-file out.parquet
```

Options:

- `-m1, --mpileup-contig-1` (required)
- `-m2, --mpileup-contig-2` (required)
- `-s, --stb-file` (required)
- `-c, --min-cov` (default: `5`)
- `-l, --min-gene-compare-len` (default: `100`)
- `-o, --output-file` (required)
- `-g, --scope` (default: `all:all`)
- `-a, --ani-method` (default: `popani`)
- `--engine` (`polars|duckdb`, default: `polars`)
- `--duckdb-memory-limit`
- `--duckdb-temp-directory`
- `--duckdb-threads`

#### `zipstrain utilities to-complete-table`

```bash
zipstrain utilities to-complete-table \
  --profile-db profile_db.parquet \
  --comp-db-file current_compare.parquet \
  --output-file remaining_pairs.csv
```

Options:

- `--profile-db` (required)
- `--comp-db-file` (optional)
- `-o, --output-file` (required)

Output columns:

- `sample_name_1`
- `sample_name_2`
- `profile_location_1`
- `profile_location_2`

Notes:

- this command does not need `--scope`, `--min-cov`, `--min-gene-compare-len`, or `--stb-file`
- it only compares the sample-pair universe implied by the profile DB against the pairs already present in the current genome comparison parquet

#### `zipstrain utilities build-genome-db`

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table sylph_abundance.tsv \
  --cache-dir genome_cache \
  --output-dir .
```

Important options:

- `--download-retries` (default: `8`)
- `--retry-backoff-seconds` (default: `10.0`)
- `--download-workers` (default: `1`)

#### Automatically build a reference database

This guide shows how to build a reference bundle directly from a Sylph abundance table using:

```bash
zipstrain utilities build-genome-db
```

Note- This functionality can also be performed at-scale with the helper program [MetaTrawl](https://github.com/OlmLab/MetaTrawl)

##### What this does

Given a Sylph abundance table, ZipStrain will:

1. Keep genomes with non-zero abundance in at least one sample.
2. Resolve/download those genomes into a persistent cache directory.
3. Reuse genomes already present in that cache (no redownload).
4. Write a concatenated reference FASTA and STB file in your output directory.

##### Command

```bash
zipstrain utilities build-genome-db \
  --tool sylph \
  --abundance-table /path/to/sylph_abundance.csv \
  --cache-dir /path/to/genome_cache \
  --output-dir /path/to/reference_bundle \
  --download-retries 8 \
  --retry-backoff-seconds 10.0 \
  --download-workers 1
```

##### Inputs

- `--tool`: abundance parser to use (currently `sylph`).
- `--abundance-table`: `.csv`, `.tsv`, or `.parquet` table.
- `--cache-dir`: persistent genome cache directory.
- `--output-dir`: output directory for the final reference bundle.

For Sylph tables, ZipStrain extracts genome accessions from the `Genome_file` column
(case-insensitive), including GTDB-style paths such as:
`gtdb_genomes_reps_r220/database/GCA/.../GCA_949068525.1_genomic.fna.gz`.

If `Genome_file` points to a local file (absolute path or path relative to the abundance-table directory),
ZipStrain loads it directly into cache first, then only downloads what is still missing.

##### Outputs

The command writes:

- `/path/to/reference_bundle/reference_genomes.fna`
- `/path/to/reference_bundle/reference_genomes.stb`
- `/path/to/reference_bundle/genome_db_build_report.txt`

##### STB format

`reference_genomes.stb` has two columns (tab-separated, no header):

- scaffold ID in the concatenated FASTA
- genome ID

Genome IDs are accessions (for example `GCF_000001405.40`).

##### Cache behavior

Inside `--cache-dir`, ZipStrain maintains:

- a local DB index (`.genome_db.parquet`)
- downloaded genomes under `genomes/`

Re-running with the same cache directory avoids redownloading genomes that already exist.

##### Retry behavior

For genomes that are not available locally/in-cache, ZipStrain retries each download with exponential backoff (default: up to 8 attempts per genome).  
If a genome still fails after retries, it is skipped, and the reference bundle is built from successfully fetched genomes.
Parallelism for remote fetch is controlled with `--download-workers`.

If you see repeated `Too Many Requests` errors on large runs:

- lower `--download-workers` (for example `1` or `2`)
- increase `--download-retries` (for example `8`)
- increase `--retry-backoff-seconds` (for example `10` to `20`)

##### Console summary

`build-genome-db` prints a short run summary:

- selected genomes (non-zero abundance)
- genomes already cached before the run
- new download attempts
- downloaded now / failed
- genomes available in cache after the run

The same summary is saved to `genome_db_build_report.txt` and includes explicit failed accession IDs (with error messages) when downloads fail.

#### `zipstrain utilities build-matrix-db`

```bash
zipstrain utilities build-matrix-db \
  --profile-dir profiles \
  --output-file matrix_db.h5 \
  --bed-file genomes_bed_file.bed \
  --stb-file reference.stb \
  --gene-range-table gene_range_table.tsv \
  --memory-limit-gb 16
```

What it does:

- scans a directory of standard ZipStrain profile parquets
- builds one matrix store directly from those profiles
- uses the BED and STB files as the explicit scaffold/genome contract for the store
- stores each genome as one sample-major dense dataset with shape `samples x positions x 4`
- positions with total coverage below `5` are zeroed during matrix build
- can optionally store scaffold-relative gene ranges for later gene ANI
- is intended for repeated cohort-scale comparison runs against the same reference set

Important options:

- `-p, --profile-dir` (required)
- `-o, --output-file` (required)
- `-g, --genome` optional genome scope (default: `all`)
- `-b, --bed-file` (required) BED file defining scaffold coordinate extents for the matrix contract
- `-s, --stb-file` (required) STB file defining scaffold-to-genome membership for the matrix contract
- `--gene-range-table` optional headerless TSV of `gene, scaffold, start, end` for gene ANI support
- `--count-dtype` stored matrix dtype (`uint16|uint32`, default: `uint16`)
- `--memory-limit-gb` approximate maximum memory budget for the entire build process (default: `16.0`)
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB (default: `128.0`)
- `--sparse` store genome matrices sparsely in HDF5

Notes:

- the output matrix store is intended for `zipstrain utilities matrix-compare`
- new matrix stores are append-friendly on the sample axis
- every input profile is interpreted against the BED+STB contract you provide here
- install matrix support with `pip install "zipstrain[matrix]"`
- the CLI shows a progress bar in an interactive terminal
- in non-interactive runs, the CLI emits throttled structured progress lines to stderr for log files
- if `--gene-range-table` is omitted, matrix compare can still compute genome ANI and IBS, but not gene ANI
- `--sparse` reduces on-disk HDF5 size, but matrix compare currently materializes sparse storage back into dense arrays when loading for comparison

#### `zipstrain utilities append-matrix-db`

```bash
zipstrain utilities append-matrix-db \
  --profile-dir new_profiles \
  --matrix-db-file matrix_db.h5 \
  --memory-limit-gb 16
```

What it does:

- scans a directory of new standard ZipStrain profile parquets
- validates that they match the existing matrix-store contract
- appends new sample rows and whole-genome matrices into the existing matrix store
- materializes newly encountered genomes when they are still compatible with the stored BED+STB contract
- ignores genomes that fall outside the stored contract and reports the ignored count

Important options:

- `-p, --profile-dir` (required)
- `-m, --matrix-db-file` (required)
- `--memory-limit-gb` approximate maximum memory budget for the append process (default: `16.0`)
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB used when rewriting an older fixed-size store (default: `128.0`)

Append behavior:

- sample names must be new
- known scaffolds and coordinate ranges must stay within the stored contract
- compatible genomes can be appended even if no matrix dataset existed for them yet
- genomes outside the stored contract are skipped and counted in the summary output

#### `zipstrain utilities matrix-db-to-hdf5`

```bash
zipstrain utilities matrix-db-to-hdf5 \
  --matrix-db-file matrix_db.duckdb \
  --output-file matrix_db.h5
```

What it does:

- converts an existing DuckDB matrix database into the current matrix-store layout
- preserves sample, genome, and scaffold metadata
- is only needed when you already have a DuckDB-based matrix database from an older workflow

Important options:

- `-m, --matrix-db-file` (required)
- `-o, --output-file` optional output HDF5 path; defaults to the same basename with `.h5`
- `--export-batch-mb` approximate matrix-store sample-axis chunk target size in MiB (default: `128.0`)

#### `zipstrain utilities matrix-compare`

```bash
zipstrain utilities matrix-compare \
  --matrix-db-file matrix_db.h5 \
  --output-file matrix_compare.duckdb \
  --memory-limit-gb 16 \
  --anchor-queue-size 1 \
  --target-queue-size 1 \
  --result-transfer-batch-size 1 \
  --loader-executor thread \
  --writer-executor thread \
  --calculate ani+ibs \
  --backend numpy
```

What it does:

- reads a per-sample genome-matrix store from `build-matrix-db` or `matrix-db-to-hdf5`
- writes results into a DuckDB compare database
- if the compare DB already exists, only pairs not yet marked completed are processed
- loads one anchor sample plus as many target samples as fit the memory budget
- computes genome ANI from dense whole-genome matrices
- computes IBS from the shared-allele boolean mask
- computes gene ANI when gene annotations are present in the matrix store and `gene` is requested
- stores result rows and completion metadata in the compare DB incrementally

Important options:

- `-m, --matrix-db-file` (required)
- `-o, --output-file` (required)
- `-g, --genome` optional genome scope (default: `all`)
- `--memory-limit-gb` approximate compare memory budget
- `--anchor-queue-size` number of torch anchor matrices to keep queued in host RAM while still transferring only one anchor at a time to the GPU (default: `1`)
- `--target-queue-size` number of torch target blocks to keep queued in host RAM; `1` preserves the current synchronous target-load behavior (default: `1`)
- `--result-transfer-batch-size` number of torch compare units to batch before transferring result vectors back to CPU (default: `1`)
- `--loader-executor` executor kind for torch loader prefetch work (`thread|process`, default: `thread`)
- `--writer-executor` executor kind for torch result writing/checkpoint work (`thread|process`, default: `thread`)
- `--calculate` matrix metrics to compute:
  - `ani`
  - `ani+ibs`
  - `+gene` or `gene` for gene ANI
  - `all` means `ani+ibs`, and also `gene` when the matrix store contains gene annotations
- `--backend` compute backend:
  - `numpy`
  - `torch`
  - `torch-cpu`
  - `torch-cuda`
  - `torch-mps`

Notes:

- install Torch support with `pip install "zipstrain[matrix]"`
- `--backend numpy` works without Torch and is the simplest CPU-only path
- on Apple Silicon, the standard `torch` wheel can use MPS
- on Linux with NVIDIA GPUs, replace Torch with the CUDA wheel that matches your system, for example:

  ```
  pip install "zipstrain[matrix]"
  pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu124
  ```

- MPS requires native macOS; Linux containers cannot expose Apple Metal
- `torch` auto-selects CUDA, then MPS, then CPU
- for torch backends, GPU work stays on the main process while loader and writer stages can run through either thread or process executors
- the compare database is resumable; rerunning the same command on the same output file only processes unfinished sample pairs
- the CLI shows a progress bar in an interactive terminal
- in non-interactive runs, the CLI emits throttled structured progress lines to stderr for log files

#### `zipstrain utilities matrix-compare-export`

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file matrix_compare.duckdb \
  --output-file matrix_compare.parquet
```

Export the gene table instead:

```bash
zipstrain utilities matrix-compare-export \
  --matrix-compare-db-file matrix_compare.duckdb \
  --output-file matrix_compare_gene.parquet \
  --table gene
```

What it does:

- reads `matrix_compare_results` from a matrix compare DuckDB
- exports the standard compare columns to parquet
- uses the stored compare metadata to choose the correct output columns
- can also export `matrix_compare_gene_results` with `--table gene`
- raises an error if `--table gene` is requested but the compare DB has no gene table

Important options:

- `-m, --matrix-compare-db-file` (required)
- `-o, --output-file` (required)
- `--table genome|gene` (optional, default `genome`)

#### `zipstrain utilities merge_parquet`

```bash
zipstrain utilities merge_parquet \
  --input-dir comps \
  --output-file merged_comparisons.parquet \
  --batch-size 5000
```

Notes:

- `--batch-size -1` keeps the current single-pass merge behavior.
- Positive `--batch-size` values first merge input files batch-by-batch into a temporary directory, then do one final lazy merge over those batch outputs.
- Progress is logged as active line-oriented batch updates and flushed immediately, which is easier to follow in cluster logs.
- When input files contain ZipStrain compare metadata, those metadata fields must match across inputs unless `--allow-mismatch` is used.
- With `--allow-mismatch`, mismatched compare metadata are rewritten to `NA` in the merged parquet metadata.

#### `zipstrain utilities build-profile-db`

```bash
zipstrain utilities build-profile-db \
  --profile-db-csv profiles.csv \
  --output-file profile_db.parquet
```

Input CSV columns:

- `profile_name`
- `profile_location`

Notes:

- By default, ZipStrain checks that all listed profiles carry matching embedded contract metadata.
- Use `--allow-mismatch` to skip that validation and build a mixed profile DB intentionally.
- The output parquet stores at least `profile_name` and `profile_location`, plus shared metadata fields derived from the listed profiles.

#### `zipstrain utilities get-coverage-stats`

```bash
zipstrain utilities get-coverage-stats \
  --profile-parquet sample_profile.parquet \
  --gene-bed reference_genomes_gene_ranges.tsv \
  --genome-bed genomes_bed_file.bed \
  --output-dir stats \
  --prefix sample1
```

What it does:

- rebuilds coverage-only genome and gene stats from an existing profile parquet
- writes:
  - `<output-dir>/<prefix>_gene_stats.parquet`
  - `<output-dir>/<prefix>_genome_stats.parquet`
- does not require read-location files
- uses the profile’s existing `gene` and `genome` columns for counts
- uses the supplied gene/genome BED files only to calculate lengths

Output columns:

- gene stats:
  - `genome`
  - `gene`
  - `length`
  - `breadth`
  - `coverage`
  - `5x_cov_sites`
  - `ber`
- genome stats:
  - `genome`
  - `length`
  - `breadth`
  - `coverage`
  - `5x_cov_sites`
  - `ber`

Important options:

- `-p, --profile-parquet` (required)
- `-g, --gene-bed` (required) — supports 4 columns `gene, scaffold, start, end` or 5 columns with genome appended
- `-b, --genome-bed` (required) — supports 3 columns `scaffold, start, end` or 4 columns with genome appended
- `-o, --output-dir` (required)
- `--prefix` (required)

#### Other Utility Commands

Use `--help` on each command for full option details:

```bash
zipstrain utilities build-null-model --help
zipstrain utilities merge_parquet --help
zipstrain utilities get-coverage-stats --help
zipstrain utilities process_mpileup --help
zipstrain utilities make_bed --help
zipstrain utilities get_genome_lengths --help
zipstrain utilities generate-genome-pairs --help
zipstrain utilities chunk-genome-compare --help
zipstrain utilities merge-stat-tables --help
zipstrain utilities strain_heterogeneity --help
zipstrain utilities build-profile-db --help
zipstrain utilities build-matrix-db --help
zipstrain utilities matrix-compare --help
zipstrain utilities presence-profile --help
zipstrain utilities process-read-locs --help
zipstrain utilities generate_stb --help
zipstrain utilities gene-range-table --help
zipstrain test --help
```


## Nextflow Implementation

This page reflects the current `zipstrain.nf` workflow in this repository.

TODO: Add a description of why you'd choose the nextflow implementation over the command line implementation

Important scope note:

- the standard profile-based compare workflows are available in Nextflow
- the newer matrix-store workflow is currently CLI-driven and is documented in [Tutorial](./Tutorial.md) and [CLI](./cli.md)
- if you want a worked standard example, the [Tutorial](./Tutorial.md) now includes both a Python/CLI route and a Nextflow route

### General Usage Instructions

TODO: Here's some instructions on how to install Nextflow, what the benefits of nextflow are, and how to set up your config file for this


#### General Running Pattern

```bash
nextflow run zipstrain.nf \
  --mode <mode> \
  --input_table <path/to/input.csv> \
  --output_dir <path/to/output_dir> \
  -c conf.config \
  -profile <docker|alpine|gutbot|blanca|fiji> \
  -resume
```

`conf.config` already defines resources for the current process set and includes example execution profiles.
The `fiji` profile is configured for Slurm plus Singularity, which is useful on clusters that provide Singularity rather than Apptainer.
Review the container paths or tags in your config before running on a new system.

#### ZipStrain Nextflow Pipeline Commands

- Read mapping with Bowtie2 (`map_reads`)
- Profile generation from BAM files (`profile`)
- End-to-end SRA to profile (`from_sra_to_profile`)
- Pairwise genome comparison across profiles (`compare_genomes`)
- Pairwise gene comparison across profiles (`compare_genes`)

#### Key Pipeline Parameters

- `--mode`: `map_reads`, `profile`, `from_sra_to_profile`, `compare_genomes`, `compare_genes`
- `--input_type`: depends on mode (`local`, `sra`, `profile_table`, `pair_table`)
- `--parallel_mode`: `single` or `batched` for comparison workflows
- `--batch_size`: number of pairs per batch when `--parallel_mode batched`
- `--batch_compare_n_parallel`: parallel jobs inside each batched comparison task
- `--compare_genome_scope`: genome scope for genome comparisons (`all` or genome ID)
- `--compare_genome_calculate`: genome metrics to compute (`all` or `ani`)
- `--compare_gene_scope`: gene scope for gene comparisons (`all:all`, `<genome>:all`, `all:<gene>`, `<genome>:<gene>`)
- `--compare_ani_method`: ANI method forwarded to compare tasks (`popani`, `conani`, `cosani_<threshold>`)
- `--compare_engine`: comparison engine for standard compare tasks (`polars` or `duckdb`). Default: `polars`
- `--compare_duckdb_memory_limit`: forwarded to single compare commands
- `--compare_calculate`: genome metrics for genome compare (`ani`, `ibs`, `identical_genes`, `all`, or `+` combinations). Default: `all`

#### Specific Notes

- For `--input_type profile_table`, use `sample_name` and `profile_location`.
- For `--input_type pair_table`, use `sample_name_1`, `sample_name_2`, `profile_location_1`, and `profile_location_2`.
- `--compare_duckdb_memory_limit` is only relevant when `--compare_engine duckdb`.
- The current Nextflow profiling workflow uses the default profiling read filters unless you edit `zipstrain.nf` directly.
- For auto-built references, genome selection comes from the merged Sylph abundance table through `zipstrain utilities build-genome-db`.


### Command 1: Map Reads (`mode=map_reads`)

#### Input Table (`--input_type local`)

Paired-end:

```csv
sample_name,reads1,reads2
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz
```

Single-end:

```csv
sample_name,reads1
S1,/data/S1.fastq.gz
S2,/data/S2.fastq.gz
```

#### Input Table (`--input_type sra`)

```csv
Run
SRR12345678
SRR12345679
```

#### Option 1.A: Use Existing Reference Genome

```bash
nextflow run zipstrain.nf \
  --mode map_reads \
  --input_type local \
  --input_table reads.csv \
  --reference_genome reference_genomes.fna \
  --stb reference_genomes.stb \
  --output_dir out_map \
  -c conf.config \
  -profile docker \
  -resume
```

Optional:

- `--index_files` to reuse existing Bowtie2 index files
- `--bowtie2_non_competitive_mapping true` to pass `-a` to Bowtie2

#### Option 1.B: Build Reference from Sylph Automatically

If `--reference_genome` is not provided, the pipeline does:

1. per-sample `sylph profile`
2. merge all per-sample Sylph abundance tables
3. `zipstrain utilities build-genome-db --tool sylph ...`
4. `prodigal` gene prediction on the generated reference FASTA
5. Bowtie2 indexing
6. mapping

```bash
nextflow run zipstrain.nf \
  --mode map_reads \
  --input_type local \
  --input_table reads.csv \
  --output_dir out_map \
  --genome_db_cache_dir genome_cache \
  --sylph_db /path/to/custom.syldb \
  -c conf.config \
  -profile docker \
  -resume
```

If `--sylph_db` is omitted, `--sylph_db_link` is used for download.

#### Map Outputs

- BAM files: `<output_dir>/*.bam`
- Sylph tables: `<output_dir>/sylph_abundance/`
- Built reference bundle (when auto-built): `<output_dir>/db_from_sylph/`

### Command 2: Generate Profiles from BAM (`mode=profile`)

#### Input Table

```csv
sample_name,bamfile
S1,/data/S1.bam
S2,/data/S2.bam
```

#### Command

```bash
nextflow run zipstrain.nf \
  --mode profile \
  --input_table bams.csv \
  --reference_genome reference_genomes.fna \
  --gene_file reference_genomes_gene.fasta \
  --stb reference_genomes.stb \
  --output_dir out_profile \
  -c conf.config \
  -profile docker \
  -resume
```

#### Profile Outputs

- `<output_dir>/*_profile.parquet`
- `<output_dir>/*_genome_stats.parquet`
- `<output_dir>/*_gene_stats.parquet`

Because the current Nextflow profiling modes pass a reference FASTA into profiling, these outputs normally include the reference-aware fields:

- profiles include `ref_base_bitmask`
- gene/genome stat tables include `ref_ani`

### Command 3: End-to-End SRA to Profile (`mode=from_sra_to_profile`)

#### Input Table

```csv
Run
SRR12345678
SRR12345679
```

#### Command

```bash
nextflow run zipstrain.nf \
  --mode from_sra_to_profile \
  --input_table sra.csv \
  --reference_genome reference_genomes.fna \
  --gene_file reference_genomes_gene.fasta \
  --stb reference_genomes.stb \
  --output_dir out_sra_profile \
  -c conf.config \
  -profile docker \
  -resume
```

#### Outputs

- `<output_dir>/profiles/*_profile.parquet`
- `<output_dir>/profiles/*_genome_stats.parquet`
- `<output_dir>/profiles/*_gene_stats.parquet`

### Command 4: Compare Genomes (`mode=compare_genomes`)

#### Input Option 4.A: All-vs-All from Profile List (`--input_type profile_table`)

```csv
sample_name,profile_location
S1,/profiles/S1_profile.parquet
S2,/profiles/S2_profile.parquet
S3,/profiles/S3_profile.parquet
```

#### Input Option 4.B: Explicit Pairs (`--input_type pair_table`)

```csv
sample_name_1,sample_name_2,profile_location_1,profile_location_2
S1,S2,/profiles/S1_profile.parquet,/profiles/S2_profile.parquet
S1,S3,/profiles/S1_profile.parquet,/profiles/S3_profile.parquet
```

#### Command

```bash
nextflow run zipstrain.nf \
  --mode compare_genomes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference_genomes.stb \
  --compare_engine polars \
  --compare_genome_scope all \
  --compare_calculate ani+ibs+identical_genes \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --output_dir out_compare_genomes \
  -c conf.config \
  -profile docker \
  -resume
```

### Command 5: Compare Genes (`mode=compare_genes`)

Input-table formats are the same as genome compare (`profile_table` or `pair_table`).

#### Command

```bash
nextflow run zipstrain.nf \
  --mode compare_genes \
  --input_type profile_table \
  --input_table profiles.csv \
  --stb reference_genomes.stb \
  --compare_engine polars \
  --compare_gene_scope all:all \
  --compare_ani_method popani \
  --parallel_mode batched \
  --batch_size 1000 \
  --batch_compare_n_parallel 4 \
  --output_dir out_compare_genes \
  -c conf.config \
  -profile docker \
  -resume
```

#### Comparison Outputs

- Final merged table: `<output_dir>/merged_comparisons.parquet`
- Intermediate batched outputs (when `parallel_mode=batched`): `<output_dir>/batch_comparisons/`
