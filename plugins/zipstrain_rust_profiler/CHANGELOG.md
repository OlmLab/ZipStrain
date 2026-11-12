# Rust profiler plugin changelog

## 0.2.0b1

- Require ZipStrain V2 and write gene-free position tables directly, including worker chunks. Remove native gene-range loading, per-position annotation, and redundant gene-statistics production. The shared V2 finalizer creates gene statistics and optional codon/selection outputs.
- Remove the native `--gene-range` option; pass gene assets to ZipStrain instead.
- Add a `codons` subcommand and `RustProfileBackend.extract_codons`. `--collect-codons` now extracts codons natively, about 28 times faster in testing and with identical counts. ZipStrain falls back to Python with older plugin builds.
- Fixed: BAQ now uses mpileup's default extended mode, and read gaps include the stretch to each scaffold's end. Rust and Python base counts, codon tables, and statistics now match. Base counts and the statistics built on them can change from earlier plugin versions wherever mates disagree or bases sit near mismatches.
- Preserve `reads_mapped` when gap uniformity is undefined, including at high coverage.

## 0.1.2

Preserve `reads_mapped` when gap uniformity is undefined, including at high coverage.

## 0.1.1

Compared with `0.1.0`, the profiler now builds a missing BAM index before starting its workers. Existing BAI and CSI indices are reused.
