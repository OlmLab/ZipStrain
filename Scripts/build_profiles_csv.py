#!/usr/bin/env python3
"""Build the profiles.csv that `zipstrain.nf --mode compare` takes as input.

Profile runs write one folder per sample containing `<sample>_profile.parquet`,
so this walks a run directory, finds those files, and emits the two-column table
the pipeline expects:

    sample_name,profiles
    SRR12324251,/abs/path/SRR12324251_profile.parquet

Paths are written absolute because Nextflow stages them from the launch
directory, which is not necessarily where this script was run.

Usage
-----
    python Scripts/build_profiles_csv.py <profile_run_dir> -o profiles.csv

    # several run directories at once
    python Scripts/build_profiles_csv.py run_a run_b -o profiles.csv

    # keep only samples listed in a file (one name per line)
    python Scripts/build_profiles_csv.py run_dir --only samples.txt -o profiles.csv
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

PROFILE_SUFFIX = "_profile.parquet"


def find_profiles(roots: list[pathlib.Path]) -> dict[str, pathlib.Path]:
    """Map sample name -> profile parquet, searching each root recursively.

    A later root wins on duplicate sample names, and the collision is reported so
    a stale run directory does not silently shadow a newer one.
    """
    found: dict[str, pathlib.Path] = {}
    for root in roots:
        if not root.exists():
            sys.exit(f"error: no such directory: {root}")
        for path in sorted(root.rglob(f"*{PROFILE_SUFFIX}")):
            sample = path.name[: -len(PROFILE_SUFFIX)]
            previous = found.get(sample)
            if previous is not None and previous != path:
                print(
                    f"warning: '{sample}' found twice; using {path} over {previous}",
                    file=sys.stderr,
                )
            found[sample] = path.resolve()
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dirs", nargs="+", type=pathlib.Path, help="Profile run directory/directories to scan.")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path("profiles.csv"), help="Output CSV (default: profiles.csv).")
    parser.add_argument("--only", type=pathlib.Path, default=None, help="File of sample names (one per line) to restrict the table to.")
    args = parser.parse_args()

    profiles = find_profiles(args.run_dirs)
    if not profiles:
        sys.exit(f"error: no *{PROFILE_SUFFIX} files found under: " + ", ".join(str(d) for d in args.run_dirs))

    if args.only is not None:
        wanted = {line.strip() for line in args.only.read_text().splitlines() if line.strip()}
        missing = sorted(wanted - profiles.keys())
        if missing:
            sys.exit("error: requested samples have no profile: " + ", ".join(missing))
        profiles = {name: path for name, path in profiles.items() if name in wanted}

    empty = [name for name, path in profiles.items() if path.stat().st_size == 0]
    if empty:
        sys.exit("error: zero-byte profiles (incomplete run?): " + ", ".join(sorted(empty)))

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_name", "profiles"])
        for sample in sorted(profiles):
            writer.writerow([sample, str(profiles[sample])])

    count = len(profiles)
    pairs = count * (count - 1) // 2
    print(f"wrote {args.output} with {count} samples -> {pairs:,} pairwise comparisons")
    if pairs:
        print(f"  at --batch_size 100 that is {-(-pairs // 100):,} batches")


if __name__ == "__main__":
    main()
