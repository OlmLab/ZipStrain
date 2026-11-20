"""The Rust codon extractor must match ZipStrain's Python extraction exactly."""

import random
import shutil
import subprocess

import polars as pl
import pytest

from zipstrain import codons, gene_info as gi
from zipstrain_rust_profiler import RustProfileBackend

pytestmark = pytest.mark.skipif(
    not shutil.which("zipstrain-rust-profiler") or not shutil.which("samtools"),
    reason="The compiled Rust profiler and samtools are required",
)
SETTINGS = [
    dict(min_baseq=13, min_mapq=0, min_read_ani=None, read_inclusion="all-mapped"),
    dict(min_baseq=13, min_mapq=0, min_read_ani=0.95, read_inclusion="paired"),
    dict(min_baseq=25, min_mapq=20, min_read_ani=0.9, read_inclusion="proper-pairs"),
    dict(min_baseq=0, min_mapq=0, min_read_ani=0.0, read_inclusion="paired"),
]


def _sequence(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def _alignment(rng, reference, start, read_length):
    """Return a CIGAR and read with matches, mismatches, indels, skips, clips, and Ns."""
    operations, sequence, position = [], [], start
    if rng.random() < .15:
        clip = rng.randint(1, 8)
        operations.append(f"{clip}{rng.choice('SH')}")
        if operations[-1].endswith("S"):
            sequence.append(_sequence(rng, clip))
    remaining, mismatches = read_length, 0
    while remaining > 0:
        length = min(remaining, rng.randint(4, 60))
        for offset in range(length):
            base = reference[position + offset - 1] if position + offset <= len(reference) else "A"
            draw = rng.random()
            if draw < .03:
                base, mismatches = rng.choice("ACGT"), mismatches + 1
            elif draw < .035:
                base = "N"
            sequence.append(base)
        operations.append(f"{length}{rng.choice('MMM=X')}")
        position += length
        remaining -= length
        if remaining > 8:
            draw = rng.random()
            if draw < .12:
                size = rng.randint(1, 4)
                operations.append(f"{size}I")
                sequence.append(_sequence(rng, size))
                remaining -= size
                mismatches += size
            elif draw < .24:
                size = rng.randint(1, 4)
                operations.append(f"{size}D")
                position += size
                mismatches += size
            elif draw < .27:
                size = rng.randint(5, 30)
                operations.append(f"{size}N")
                position += size
    if rng.random() < .1:
        clip = rng.randint(1, 8)
        operations.append(f"{clip}S")
        sequence.append(_sequence(rng, clip))
    return "".join(operations), "".join(sequence), mismatches


def _write_bam(root, references, rng, fragments):
    rows = []
    for _ in range(fragments):
        name = rng.choice(list(references))
        reference = references[name]
        first = rng.randint(1, len(reference) - 150)
        # Reuse some names, so mates, groups, and later fragments can collide in the cache.
        qname = f"f{rng.randint(0, fragments // 2)}"
        group = rng.choice(["", "\tRG:Z:g1", "\tRG:Z:g2"])
        second = max(first, min(len(reference) - 100, first + rng.randint(-20, 90)))
        for mate, start in enumerate((first, second)):
            if mate == 1 and rng.random() < .2:
                continue
            flag = 1 | (64 if mate == 0 else 128)
            if rng.random() < .8:
                flag |= 2
            if rng.random() < .05:
                flag |= 8
            for bit in (0x100, 0x200, 0x400, 0x800):
                if rng.random() < .02:
                    flag |= bit
            cigar, sequence, mismatches = _alignment(rng, reference, start, rng.randint(30, 150))
            quality = "*" if rng.random() < .02 else "".join(
                chr(33 + rng.choice((2, 13, 20, 25, 30, 40))) for _ in sequence)
            nm = "" if rng.random() < .1 else f"\tNM:i:{mismatches}"
            rnext = "=" if rng.random() < .9 else rng.choice(list(references))
            mapq = rng.choice((0, 10, 30, 60))
            rows.append((name, start, f"{qname}\t{flag}\t{name}\t{start}\t{mapq}\t{cigar}\t{rnext}\t"
                                      f"{second if mate == 0 else first}\t0\t{sequence}\t{quality}{nm}{group}"))
    header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{n}\tLN:{len(s)}\n" for n, s in references.items())
    header += "@RG\tID:g1\n@RG\tID:g2\n"
    order = {name: index for index, name in enumerate(references)}
    rows.sort(key=lambda row: (order[row[0]], row[1]))
    sam = root/"reads.sam"
    sam.write_text(header + "\n".join(row[2] for row in rows) + "\n")
    bam = root/"reads.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    return bam


def _catalog(root, references, rng):
    """Write forward, reverse, nested, and overlapping genes."""
    records = []
    for name, sequence in references.items():
        for index in range(max(1, len(sequence) // 200)):
            start = rng.randint(1, len(sequence) - 40)
            end = min(len(sequence), start + 3 * rng.randint(4, 80) - 1)
            end -= (end - start + 1) % 3
            strand = rng.choice((1, -1))
            segment = sequence[start - 1:end]
            if strand == -1:
                segment = segment[::-1].translate(str.maketrans("ACGT", "TGCA"))
            records.append(f">{name}_{index + 1} # {start} # {end} # {strand} # partial=00\n{segment}\n")
    (root/"ref.fa").write_text("".join(f">{n}\n{s}\n" for n, s in references.items()))
    (root/"genes.fa").write_text("".join(records))
    (root/"ref.stb").write_text("".join(f"{n}\tgenome_{n}\n" for n in references))
    return gi.write_gene_info(root/"ref.fa", root/"genes.fa", root/"ref.stb", root/"gene_info.parquet")


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("settings", SETTINGS, ids=lambda value: f"{value['read_inclusion']}-{value['min_baseq']}")
def test_rust_extraction_matches_python(tmp_path, seed, settings):
    rng = random.Random(seed)
    references = {name: _sequence(rng, length) for name, length in (("s1", 3000), ("s2", 1200), ("nogenes", 400))}
    catalog = _catalog(tmp_path, {n: s for n, s in references.items() if n != "nogenes"}, rng)
    bam = _write_bam(tmp_path, references, rng, 1500)
    genes = gi.read_gene_info(catalog)
    python, rust = tmp_path/"python.parquet", tmp_path/"rust.parquet"
    codons._extract_raw_codons(bam, genes, python, **settings)
    RustProfileBackend().extract_codons(bam, genes, rust, threads=2, **settings)
    keys = ["gene_id", "codon_index", "codon"]
    expected = pl.read_parquet(python).sort(keys)
    actual = pl.read_parquet(rust).cast(codons.COUNT_SCHEMA).sort(keys)
    assert expected.height > 50
    assert actual.equals(expected)


def test_rust_extraction_rejects_unsorted_bam(tmp_path):
    rng = random.Random(9)
    references = {"s1": _sequence(rng, 2000), "s2": _sequence(rng, 2000)}
    catalog = _catalog(tmp_path, references, rng)
    header = "".join(f"@SQ\tSN:{n}\tLN:2000\n" for n in references)

    def read(name, start):
        return f"r\t0\t{name}\t{start}\t60\t30M\t*\t0\t0\t{references[name][start - 1:start + 29]}\t{'I' * 30}\tNM:i:0"

    for label, rows in (("position", [read("s1", 100), read("s1", 50)]),
                        ("scaffold", [read("s1", 100), read("s2", 100), read("s1", 200)])):
        sam = tmp_path/f"{label}.sam"
        sam.write_text(header + "\n".join(rows) + "\n")
        bam = tmp_path/f"{label}.bam"
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
        settings = dict(min_baseq=13, min_mapq=0, min_read_ani=None, read_inclusion="all-mapped")
        with pytest.raises(ValueError, match="coordinate-sorted"):
            codons._extract_raw_codons(bam, gi.read_gene_info(catalog), tmp_path/"p.parquet", **settings)
        with pytest.raises(subprocess.CalledProcessError):
            RustProfileBackend().extract_codons(bam, gi.read_gene_info(catalog), tmp_path/"r.parquet", threads=1, **settings)
