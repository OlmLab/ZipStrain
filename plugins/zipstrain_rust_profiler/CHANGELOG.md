# Rust profiler plugin changelog

## 0.2.0.dev0 (unreleased)

- Require ZipStrain V2 and write gene-free position tables directly, including worker chunks. Remove native gene-range loading, per-position annotation, and redundant gene-statistics production. The shared V2 finalizer creates gene statistics and optional codon/selection outputs.
- Remove the native `--gene-range` option; pass gene assets to ZipStrain instead.

## 0.1.1

Compared with `0.1.0`, the profiler now builds a missing BAM index before starting its workers. Existing BAI and CSI indices are reused.
