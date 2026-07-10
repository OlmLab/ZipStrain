"""zipstrain.mapping
========================
Read-to-BAM mapping for ZipStrain.

This module powers ``zipstrain map``: it turns a table of sequencing reads into
sorted, indexed BAM files ready for ``zipstrain profile``. When no reference is
supplied it uses Sylph to pick reference genomes automatically (downloading and
caching them), builds a concatenated reference FASTA + STB, and then maps the
reads against it with Bowtie2.

The heavy lifting (Sylph, Bowtie2, samtools, prodigal) is done by shelling out
to those tools, so they must be installed and on ``PATH``.
"""
from __future__ import annotations

import csv
import os
import pathlib
import shlex
import shutil
import subprocess
import urllib.request
from dataclasses import dataclass
from typing import Callable

import polars as pl

from zipstrain import build_db as bdb


# A callback invoked with a short human-readable message as each step starts.
ProgressCallback = Callable[[str], None]


def _notify(callback: ProgressCallback | None, message: str) -> None:
    if callback is not None:
        callback(message)


# The GTDB Sylph database used by the bundled Nextflow pipeline.
DEFAULT_SYLPH_DB_URL = "http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb"

REFERENCE_FASTA_NAME = "reference_genomes.fna"
REFERENCE_STB_NAME = "reference_genomes.stb"
REFERENCE_GENE_FASTA_NAME = "reference_genomes_gene.fasta"
REFERENCE_TAXONOMY_NAME = "reference_genomes_taxonomy.tsv"
SAMPLES_TABLE_NAME = "samples.txt"
INTERMEDIATE_DIRNAME = "intermediate_files"

# GTDB r220 taxonomy tables (matching the default Sylph database). Each row maps
# an accession (RS_GCF.../GB_GCA...) to its full GTDB lineage.
GTDB_TAXONOMY_URLS = (
    "https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_taxonomy_r220.tsv.gz",
    "https://data.gtdb.ecogenomic.org/releases/release220/220.0/ar53_taxonomy_r220.tsv.gz",
)


@dataclass
class ReadSample:
    """A single sample's reads. ``reads2`` is ``None`` for single-end data."""

    sample_name: str
    reads1: pathlib.Path
    reads2: pathlib.Path | None = None

    @property
    def is_paired(self) -> bool:
        return self.reads2 is not None


def require_tool(name: str) -> str:
    """Return the resolved path to ``name`` or raise a clear error if missing."""
    resolved = shutil.which(name)
    if resolved is None:
        raise RuntimeError(
            f"Required tool '{name}' was not found on PATH. Install it (e.g. via "
            f"conda: `conda install -c bioconda {name}`) and try again."
        )
    return resolved


def _run(command: list[str], *, stdout=subprocess.PIPE, cwd: pathlib.Path | None = None) -> None:
    """Run ``command``, raising a helpful error on failure.

    By default the command's stdout is captured (and discarded on success) so
    chatty tools like ``bowtie2-build`` do not clutter the terminal; pass an open
    file handle as ``stdout`` to redirect it (e.g. for ``sylph profile``).
    """
    result = subprocess.run(command, stdout=stdout, stderr=subprocess.PIPE, cwd=cwd)
    if result.returncode != 0:
        captured_out = result.stdout.decode(errors="replace") if result.stdout else ""
        stderr = result.stderr.decode(errors="replace") if result.stderr else ""
        raise RuntimeError(
            f"Command failed ({' '.join(shlex.quote(c) for c in command)}):\n"
            f"{stderr}{captured_out}"
        )


def _run_shell(pipeline: str, *, cwd: pathlib.Path | None = None) -> None:
    """Run a shell pipeline string (used for the bowtie2 | samtools pipe)."""
    result = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        stderr=subprocess.PIPE,
        cwd=cwd,
    )
    if result.returncode != 0:
        stderr = result.stderr.decode(errors="replace") if result.stderr else ""
        raise RuntimeError(f"Command failed:\n{pipeline}\n{stderr}")


def read_reads_table(path: str | pathlib.Path) -> list[ReadSample]:
    """Parse a reads table CSV with columns ``sample_name,reads1[,reads2]``."""
    path = pathlib.Path(path)
    samples: list[ReadSample] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "sample_name" not in reader.fieldnames or "reads1" not in reader.fieldnames:
            raise ValueError(
                "Reads table must be a CSV with at least 'sample_name' and 'reads1' columns."
            )
        for row in reader:
            sample_name = (row.get("sample_name") or "").strip()
            reads1 = (row.get("reads1") or "").strip()
            reads2 = (row.get("reads2") or "").strip()
            if not sample_name or not reads1:
                raise ValueError(f"Every row needs a sample_name and reads1 (got: {row}).")
            samples.append(
                ReadSample(
                    sample_name=sample_name,
                    reads1=pathlib.Path(reads1),
                    reads2=pathlib.Path(reads2) if reads2 else None,
                )
            )
    if not samples:
        raise ValueError(f"No samples found in reads table: {path}")
    return samples


def ensure_sylph_db(sylph_db: str | pathlib.Path, url: str = DEFAULT_SYLPH_DB_URL) -> pathlib.Path:
    """Return the Sylph database path, downloading it from ``url`` if absent.

    Note: the default GTDB database is very large (~14 GB); the first download
    can take a while.
    """
    db_path = pathlib.Path(sylph_db)
    if db_path.exists() and db_path.stat().st_size > 0:
        return db_path
    db_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Sylph database not found at {db_path}; downloading from {url} (this may be large)...")
    tmp_path = db_path.with_suffix(db_path.suffix + ".part")
    with urllib.request.urlopen(url) as response, open(tmp_path, "wb") as out:
        shutil.copyfileobj(response, out)
    tmp_path.replace(db_path)
    return db_path


def run_sylph_profile(
    *, sylph_db: pathlib.Path, sample: ReadSample, output_tsv: pathlib.Path, threads: int
) -> pathlib.Path:
    """Run ``sylph profile`` for one sample, writing its abundance TSV."""
    sylph = require_tool("sylph")
    command = [sylph, "profile", str(sylph_db), "-t", str(threads)]
    if sample.is_paired:
        command += ["-1", str(sample.reads1), "-2", str(sample.reads2)]
    else:
        command += ["-U", str(sample.reads1)]
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "wb") as handle:
        _run(command, stdout=handle)
    return output_tsv


def merge_sylph_abundances(tsvs: list[pathlib.Path], output_tsv: pathlib.Path) -> pathlib.Path:
    """Concatenate per-sample Sylph abundance TSVs, keeping a single header."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w") as out:
        header_written = False
        for tsv in tsvs:
            with open(tsv) as handle:
                lines = handle.readlines()
            if not lines:
                continue
            if header_written:
                lines = lines[1:]  # skip repeated header
            else:
                header_written = True
            out.writelines(lines)
    return output_tsv


def predict_genes(*, reference_fasta: pathlib.Path, output_gene_fasta: pathlib.Path) -> pathlib.Path:
    """Predict genes on the reference FASTA with prodigal (meta mode)."""
    prodigal = require_tool("prodigal")
    _run([prodigal, "-i", str(reference_fasta), "-d", str(output_gene_fasta), "-p", "meta"])
    return output_gene_fasta


def build_bowtie2_index(*, reference_fasta: pathlib.Path, index_prefix: pathlib.Path, threads: int) -> pathlib.Path:
    """Build a Bowtie2 index for ``reference_fasta`` at ``index_prefix``."""
    bowtie2_build = require_tool("bowtie2-build")
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    _run([bowtie2_build, "--threads", str(threads), str(reference_fasta), str(index_prefix)])
    return index_prefix


def map_sample(
    *,
    sample: ReadSample,
    index_prefix: pathlib.Path,
    output_bam: pathlib.Path,
    threads: int,
    non_competitive: bool = False,
) -> pathlib.Path:
    """Map one sample's reads with Bowtie2, producing a sorted, indexed BAM."""
    bowtie2 = require_tool("bowtie2")
    samtools = require_tool("samtools")
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    q = shlex.quote
    competitiveness = "-a" if non_competitive else ""
    if sample.is_paired:
        reads_args = f"-1 {q(str(sample.reads1))} -2 {q(str(sample.reads2))}"
    else:
        reads_args = f"-U {q(str(sample.reads1))}"

    pipeline = (
        f"{q(bowtie2)} -x {q(str(index_prefix))} {reads_args} {competitiveness} "
        f"--threads {threads} "
        f"| {q(samtools)} view -bS -F 4 - "
        f"| {q(samtools)} sort -@ {threads} -o {q(str(output_bam))} -"
    )
    _run_shell(pipeline)
    _run([samtools, "index", str(output_bam)])
    return output_bam


def build_reference_with_sylph(
    *,
    samples: list[ReadSample],
    sylph_db: pathlib.Path,
    genome_cache_dir: pathlib.Path,
    output_dir: pathlib.Path,
    intermediate_dir: pathlib.Path,
    threads: int,
    progress_callback: ProgressCallback | None = None,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Run Sylph per sample, merge, and build a reference FASTA + STB.

    Returns ``(reference_fasta, stb_file)`` written into ``output_dir``.
    """
    abundance_dir = intermediate_dir / "sylph_abundance"
    per_sample_tsvs = []
    for index, sample in enumerate(samples, start=1):
        _notify(progress_callback, f"Sylph profiling {sample.sample_name} ({index}/{len(samples)})")
        tsv = run_sylph_profile(
            sylph_db=sylph_db,
            sample=sample,
            output_tsv=abundance_dir / f"{sample.sample_name}_sylph_abundance.tsv",
            threads=threads,
        )
        per_sample_tsvs.append(tsv)

    _notify(progress_callback, "Building reference from Sylph abundances (downloading genomes as needed)")
    merged = merge_sylph_abundances(per_sample_tsvs, abundance_dir / "sylph_abundance.tsv")

    out_fasta, out_stb, _extracted, _report, _summary = bdb.build_reference_from_abundance(
        tool_name="sylph",
        abundance_table=merged,
        cache_dir=genome_cache_dir,
        output_dir=output_dir,
        download_workers=max(1, threads),
    )
    out_fasta = pathlib.Path(out_fasta)
    if not out_fasta.exists() or out_fasta.stat().st_size == 0:
        raise RuntimeError(
            "Sylph did not detect any reference genomes in these reads, so the "
            "reference is empty and there is nothing to map against. This usually "
            "means the reads are too shallow, or the organisms are absent from the "
            f"Sylph database. Per-sample abundance tables are in {abundance_dir}. "
            "Try deeper reads, a different Sylph database, or provide "
            "--reference-fasta / --stb-file directly."
        )

    # Best-effort GTDB taxonomy for the detected genomes. Never fail the map over it.
    try:
        _notify(progress_callback, "Annotating reference genomes with GTDB taxonomy")
        write_reference_taxonomy(
            stb_file=pathlib.Path(out_stb),
            cache_dir=pathlib.Path(genome_cache_dir),
            output_file=pathlib.Path(output_dir) / REFERENCE_TAXONOMY_NAME,
        )
    except Exception as exc:  # noqa: BLE001 - taxonomy is optional
        print(f"Warning: could not build GTDB taxonomy table ({exc}); continuing without it.")

    return out_fasta, pathlib.Path(out_stb)


def ensure_gtdb_taxonomy(cache_dir: pathlib.Path, urls=GTDB_TAXONOMY_URLS) -> list[pathlib.Path]:
    """Download the GTDB taxonomy TSVs into ``cache_dir`` (cached across runs)."""
    cache_dir = pathlib.Path(cache_dir)
    tax_dir = cache_dir / "gtdb_taxonomy"
    tax_dir.mkdir(parents=True, exist_ok=True)
    local_paths = []
    for url in urls:
        dest = tax_dir / url.rsplit("/", 1)[-1]
        if not (dest.exists() and dest.stat().st_size > 0):
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url) as response, open(tmp, "wb") as out:
                shutil.copyfileobj(response, out)
            tmp.replace(dest)
        local_paths.append(dest)
    return local_paths


def write_reference_taxonomy(
    *, stb_file: pathlib.Path, cache_dir: pathlib.Path, output_file: pathlib.Path
) -> pathlib.Path | None:
    """Write a ``genome\\tgenome_taxonomy`` table for the reference's genomes.

    Genome IDs are the accessions in the STB; they are matched to GTDB lineages
    (dropping the ``RS_``/``GB_`` prefix). Returns the output path, or ``None``
    when no genome could be matched.
    """
    genome_ids = (
        pl.scan_csv(stb_file, separator="\t", has_header=False)
        .select(pl.col("column_2").alias("genome"))
        .unique()
        .collect()
        .get_column("genome")
        .to_list()
    )
    wanted = {str(g).strip() for g in genome_ids}

    tax_paths = ensure_gtdb_taxonomy(cache_dir)
    lineage: dict[str, str] = {}
    for tax_path in tax_paths:
        table = pl.read_csv(
            tax_path, separator="\t", has_header=False, new_columns=["gtdb_accession", "taxonomy"]
        ).with_columns(
            pl.col("gtdb_accession").str.replace(r"^(RS_|GB_)", "").alias("genome")
        )
        for genome, taxonomy in zip(table.get_column("genome"), table.get_column("taxonomy")):
            if genome in wanted:
                lineage[genome] = taxonomy

    if not lineage:
        return None
    output_file.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(lineage.items())
    pl.DataFrame({"genome": [g for g, _ in rows], "genome_taxonomy": [t for _, t in rows]}).write_csv(
        output_file, separator="\t"
    )
    return output_file


def run_map(
    *,
    reads_table: str | pathlib.Path,
    output_dir: str | pathlib.Path,
    reference_fasta: str | pathlib.Path | None = None,
    stb_file: str | pathlib.Path | None = None,
    sylph_db: str | pathlib.Path | None = None,
    sylph_db_url: str = DEFAULT_SYLPH_DB_URL,
    genome_cache_dir: str | pathlib.Path | None = None,
    threads: int = 4,
    predict_genes_flag: bool = False,
    non_competitive: bool = False,
    progress_callback: ProgressCallback | None = None,
) -> dict[str, pathlib.Path]:
    """End-to-end map: reads -> reference (Sylph or supplied) -> sorted BAMs.

    Returns a dict with keys ``reference_fasta``, ``stb_file``, ``samples_table``
    (and ``gene_fasta`` when ``predict_genes_flag`` is set). ``progress_callback``
    is invoked with a short message as each step begins.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    intermediate_dir = output_dir / INTERMEDIATE_DIRNAME
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    samples = read_reads_table(reads_table)
    _notify(progress_callback, f"Read {len(samples)} sample(s) from the reads table")

    # 1. Resolve the reference: use the supplied one, or build it with Sylph.
    if reference_fasta is not None:
        if stb_file is None:
            raise ValueError("--stb-file is required when --reference-fasta is provided.")
        ref_fasta = pathlib.Path(reference_fasta)
        ref_stb = pathlib.Path(stb_file)
    else:
        if genome_cache_dir is None:
            raise ValueError(
                "--genome-cache-dir is required when no --reference-fasta is provided "
                "(it caches genomes downloaded during Sylph-based reference building)."
            )
        if sylph_db is None:
            raise ValueError(
                "--sylph-db is required when no --reference-fasta is provided "
                "(path to the Sylph database; it is downloaded there if missing)."
            )
        resolved_db = ensure_sylph_db(sylph_db, url=sylph_db_url)
        ref_fasta, ref_stb = build_reference_with_sylph(
            samples=samples,
            sylph_db=resolved_db,
            genome_cache_dir=pathlib.Path(genome_cache_dir),
            output_dir=output_dir,
            intermediate_dir=intermediate_dir,
            threads=threads,
            progress_callback=progress_callback,
        )

    results: dict[str, pathlib.Path] = {"reference_fasta": ref_fasta, "stb_file": ref_stb}

    # 2. Optional gene prediction.
    if predict_genes_flag:
        _notify(progress_callback, "Predicting genes with prodigal")
        results["gene_fasta"] = predict_genes(
            reference_fasta=ref_fasta,
            output_gene_fasta=output_dir / REFERENCE_GENE_FASTA_NAME,
        )

    # 3. Build the Bowtie2 index (kept in intermediate_files/bt2, not the CWD).
    _notify(progress_callback, "Building Bowtie2 index")
    index_prefix = intermediate_dir / "bt2" / ref_fasta.name
    build_bowtie2_index(reference_fasta=ref_fasta, index_prefix=index_prefix, threads=threads)

    # 4. Map every sample to a sorted, indexed BAM.
    sample_bams: list[tuple[str, pathlib.Path]] = []
    for index, sample in enumerate(samples, start=1):
        _notify(progress_callback, f"Mapping {sample.sample_name} ({index}/{len(samples)})")
        bam = map_sample(
            sample=sample,
            index_prefix=index_prefix,
            output_bam=output_dir / f"{sample.sample_name}.bam",
            threads=threads,
            non_competitive=non_competitive,
        )
        sample_bams.append((sample.sample_name, bam))

    # 5. Write a samples table ready for `zipstrain profile`.
    samples_table = output_dir / SAMPLES_TABLE_NAME
    with open(samples_table, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_name", "bamfile"])
        for sample_name, bam in sample_bams:
            writer.writerow([sample_name, str(bam)])
    results["samples_table"] = samples_table

    return results
