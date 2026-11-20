"""Build a pair of closely related strains with known (planted) mosaic blocks.

Strain A is the reference genome itself. Strain B is a copy of A that differs at a
uniform background rate (0.1% by default, about 99.9% ANI), except inside planted blocks:

* low-divergence blocks are identical to A, as after a recent recombination
  between the two lineages;
* high-divergence blocks differ at 2%, as after importing DNA from a more
  distant lineage.

Outputs (all under --outdir):
  strain_A.fna, strain_B.fna     genomes that reads are simulated from
  reference.stb                  scaffold-to-genome table for ZipStrain
  truth/planted_blocks.tsv       one row per planted block (1-based, inclusive)
  truth/snps.tsv                 every position where B differs from A
"""

import argparse
from pathlib import Path

import numpy as np

# name, kind, start, end (1-based, inclusive), substitution rate inside the block
PLANTED_BLOCKS = [
    ("R1", "low_divergence", 300_001, 350_000, 0.0),
    ("I1", "high_divergence", 700_001, 730_000, 0.02),
    ("R2", "low_divergence", 1_100_001, 1_120_000, 0.0),
    ("I2", "high_divergence", 1_500_001, 1_505_000, 0.02),
    ("R3", "low_divergence", 1_850_001, 1_855_000, 0.0),
]


def read_fasta(path):
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    raise SystemExit("This tutorial expects a single-sequence FASTA")
                header = line[1:].split()[0]
            else:
                chunks.append(line.strip().upper())
    return header, "".join(chunks)


def write_fasta(path, name, sequence, width=80):
    with open(path, "w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), width):
            handle.write(sequence[start:start + width] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, help="Single-sequence reference FASTA (strain A)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--genome-name", default="TIGR4", help="Genome name written to the STB file")
    parser.add_argument("--background-rate", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    chrom, sequence = read_fasta(args.reference)
    length = len(sequence)

    # Per-position substitution rate: background everywhere, overridden inside blocks.
    rate = np.full(length, args.background_rate)
    for _, _, start, end, block_rate in PLANTED_BLOCKS:
        rate[start - 1:end] = block_rate
    mutated = np.flatnonzero(rng.random(length) < rate)

    bases = np.frombuffer(sequence.encode(), dtype="S1").copy()
    alternatives = {b"A": b"CGT", b"C": b"AGT", b"G": b"ACT", b"T": b"ACG"}
    snps = []
    for index in mutated:
        ref = bases[index]
        alt = bytes([rng.choice(list(alternatives[ref]))])
        bases[index] = alt
        snps.append((index + 1, ref.decode(), alt.decode()))
    strain_b = bases.tobytes().decode()

    outdir = Path(args.outdir)
    (outdir / "truth").mkdir(parents=True, exist_ok=True)
    write_fasta(outdir / "strain_A.fna", chrom, sequence)
    write_fasta(outdir / "strain_B.fna", chrom, strain_b)
    (outdir / "reference.stb").write_text(f"{chrom}\t{args.genome_name}\n")

    positions = np.array([pos for pos, _, _ in snps])
    with open(outdir / "truth" / "snps.tsv", "w") as handle:
        handle.write("chrom\tpos\tref\talt\n")
        handle.writelines(f"{chrom}\t{pos}\t{ref}\t{alt}\n" for pos, ref, alt in snps)
    with open(outdir / "truth" / "planted_blocks.tsv", "w") as handle:
        handle.write("block\tkind\tchrom\tstart\tend\tlength\tsubstitution_rate\tsnps\n")
        for name, kind, start, end, block_rate in PLANTED_BLOCKS:
            count = int(((positions >= start) & (positions <= end)).sum())
            handle.write(f"{name}\t{kind}\t{chrom}\t{start}\t{end}\t{end - start + 1}\t{block_rate}\t{count}\n")

    print(f"{chrom}: {length:,} bp; strain B differs at {len(snps):,} positions "
          f"({len(snps) / length:.4%}); {len(PLANTED_BLOCKS)} planted blocks")


if __name__ == "__main__":
    main()
