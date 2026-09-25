# zipstrain_rust_profiler

A ZipStrain plugin that provides a native, multithreaded profiling backend. It
can profile up to ~10x faster than the built-in Python profiler, with better
memory management on large BAM files, while producing the same outputs. This
package requires ZipStrain 1.2.0 or newer and is selected with `--backend rust_profiler`. Without that
flag, ZipStrain continues to use its built-in Python profiler.

```bash
python -m pip install zipstrain_rust_profiler
zipstrain profile --input-table mapped/samples.txt --reference-fasta ref.fna \
  --stb-file ref.stb --run-dir profiled --backend rust_profiler
```

See the [plugin guide](https://olmlab.github.io/ZipStrain/plugins/) for inputs,
output files, filters, and Nextflow usage. Wheels contain a native executable;
source installs require a Rust toolchain and native build prerequisites.
