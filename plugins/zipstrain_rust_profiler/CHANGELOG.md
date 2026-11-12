# Rust profiler plugin changelog

## Unreleased

Preserve `reads_mapped` when gap uniformity is undefined, including at high coverage.

## 0.1.1

Compared with `0.1.0`, the profiler now builds a missing BAM index before starting its workers. Existing BAI and CSI indices are reused.
