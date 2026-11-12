# zipstrain_rust_profiler

The Rust profiler builds a missing BAM index beside the input before profiling. The BAM must be coordinate-sorted and its directory writable unless an index already exists.

A ZipStrain plugin that provides a native, multithreaded profiling backend. It
can profile up to ~10x faster than the built-in Python profiler, with better
memory management on large BAM files, while producing the same outputs. This
beta requires ZipStrain 2.0.0b1 or a later 2.x release and is selected with `--backend rust_profiler`. Without that
flag, ZipStrain continues to use its built-in Python profiler.

```bash
python -m pip install "zipstrain-rust-profiler==0.2.0b1"
zipstrain profile --input-table mapped/samples.txt --reference-fasta ref.fna \
  --stb-file ref.stb --run-dir profiled --backend rust_profiler
```

See the [plugin guide](https://olmlab.github.io/ZipStrain/plugins/) for inputs,
output files, filters, and Nextflow usage. Wheels contain a native executable;
source installs require a Rust toolchain and native build prerequisites.

For V2, Rust writes gene-free position Parquets directly: `chrom`, `genome`,
`pos`, `A`, `C`, `G`, `T`, plus `ref_base_bitmask` when a reference is provided.
Worker chunks use the same schema. Rust does not load gene annotations, assign
per-position gene labels, or write a gene-statistics file. The native executable
no longer accepts `--gene-range`.

ZipStrain's shared V2 finalizer consumes `gene_info.parquet` and creates the
15-column gene-statistics table after native profiling. With `--collect-codons`,
the plugin makes an additional BAM pass in Rust for observed codons, with counts
identical to Python extraction, and ZipStrain computes the six selection metrics
using the requested `--allele-integration` mode. Base counts, genome statistics,
and codon extraction run in Rust; gene statistics and selection metrics run in
Python/Polars.
See [the V2 method documentation](../../docs/v2.md) for formulas and exclusions.
