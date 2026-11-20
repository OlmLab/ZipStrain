# Rust profiler plugin changelog

## 0.2.0

Stable companion for ZipStrain 2.0.0, requiring core `>=2.0.0,<3`.
Profiling behavior is unchanged from 0.2.0b2.

- Write gene-free position Parquets and extract read-phased codon evidence natively; the shared V2 finalizer produces gene statistics and selection metrics.
- Match Python for BAQ, overlapping mates, trailing scaffold gaps, and coverage/breadth/BER before allele filtering. Preserve mapped-read counts when gap uniformity is undefined.
- Publish Linux x86_64 and macOS x86_64/arm64 wheels plus a source distribution. Re-profile older BAMs for corrected low-coverage statistics. Keep plugin 0.1.2 for ZipStrain 1.x.

## 0.2.0b2

- Fixed: genome `coverage`, `coverage_median`, `coverage_std`, `breadth` and `ber` are counted from the depth before the error filter, matching ZipStrain. The filter removes a position's only allele when too few reads support it, so the old values under-reported genomes below about 2x coverage.

- Require ZipStrain `>=2.0.0b2,<3` to pair the corrected native statistics with the matching core beta.
- Keep zero-count position rows so the comparison genome inventory is preserved.
- Added low-coverage Python/Rust parity coverage.

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
