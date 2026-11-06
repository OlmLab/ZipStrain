"""V2 reference annotations, identity checks, and overlap-safe interval summaries.

Coordinates are one-based and inclusive. Gene IDs are content-derived UInt64
values, never row numbers; the complete mapping is authenticated by metadata.
See ``docs/v2.md`` for the file contracts and examples of every public function.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import polars as pl


GENE_SCHEMA = {
    "gene_id": pl.UInt64, "gene": pl.String, "genome": pl.String,
    "chrom": pl.String, "start": pl.Int64, "end": pl.Int64,
    "length": pl.Int64, "strand": pl.Int8, "phase": pl.UInt8,
    "genetic_code": pl.UInt8, "partial_5p": pl.Boolean, "partial_3p": pl.Boolean,
}
CODON_SCHEMA = {"gene_id": pl.UInt64, "codon_index": pl.UInt32, "reference_codon": pl.String}
SELECTION_COLUMNS = ("dN", "dS", "dN_dS", "pN", "pS", "pN_pS")
GENE_STATS_COLUMNS = ("gene_id", "genome", "chrom", "start", "end", "length",
                      "coverage", "breadth", "ref_ani", *SELECTION_COLUMNS)
REFERENCE_ID = "zipstrain_reference_id"
ANNOTATION_ID = "zipstrain_annotation_id"
SCHEMA_KEY = "zipstrain_profile_schema"


def fasta_records(path: str | Path):
    """Yield ``(full_header, uppercase_sequence)`` from a plain FASTA.

    Only one sequence is buffered. Empty sequence records are rejected, as are
    bases before the first header. Sequence alphabet validation is left to users.
    """
    header, sequence = None, []
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not sequence:
                        raise ValueError(f"Empty FASTA sequence: {header}")
                    yield header, "".join(sequence).upper()
                header, sequence = line[1:], []
            elif header is None:
                raise ValueError("FASTA sequence before first header")
            else:
                sequence.append(line)
        if header is not None:
            if not sequence:
                raise ValueError(f"Empty FASTA sequence: {header}")
            yield header, "".join(sequence).upper()


def content_digest(value) -> str:
    """SHA256 canonical JSON: mapping keys are sorted; callers must order sequence rows."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def reference_identity(reference_fasta, stb: pl.LazyFrame) -> str:
    """Fingerprint reference sequences and their scaffold/genome assignments.

    Ignore FASTA record order, wrapping, case, and descriptions. Require exactly
    one STB assignment for every reference scaffold and reject duplicate names.
    Used to ensure explicit profiling inputs really match a supplied catalog.
    """
    frame = stb.select("scaffold", "genome").collect()
    if frame["scaffold"].n_unique() != frame.height:
        raise ValueError("Duplicate scaffold in STB")
    mapping, records, seen = dict(frame.iter_rows()), [], set()
    for header, sequence in fasta_records(reference_fasta):
        chrom = header.split()[0]
        if chrom in seen or chrom not in mapping:
            raise ValueError("Duplicate reference scaffold or missing STB entry")
        seen.add(chrom)
        records.append((chrom, mapping[chrom], len(sequence), hashlib.sha256(sequence.encode()).hexdigest()))
    if seen != set(mapping):
        raise ValueError("STB contains scaffolds absent from reference")
    return content_digest(sorted(records))


def write_gene_info(reference_fasta, gene_fasta, stb_file, output_file) -> Path:
    """Write the V2 gene catalog and adjacent ``reference_codons.parquet``.

    Parse Prodigal nucleotide FASTA headers; validate coordinates, strand, phase,
    translation table, and unique gene names. IDs include the genome, scaffold
    sequence identity, name, coordinates, strand, phase, and genetic code. A
    detected ID collision fails rather than merging genes. Tables 1 and 11 are
    currently supported. Reference codons follow the coding orientation and omit
    incomplete terminal codons. With no gene FASTA, write typed empty catalogs.
    """
    from zipstrain.profile import read_stb

    mapping = read_stb(stb_file).collect()
    if mapping["scaffold"].n_unique() != mapping.height:
        raise ValueError("Duplicate scaffold in STB")
    genomes = dict(mapping.iter_rows())
    annotations = {}
    names = set()
    if gene_fasta is not None:
        for header, _sequence in fasta_records(gene_fasta):
            parts = [part.strip() for part in header.split("#")]
            if len(parts) < 4:
                raise ValueError("Gene FASTA must contain Prodigal coordinate headers")
            name = parts[0].split()[0]
            if name in names:
                raise ValueError(f"Duplicate gene name: {name}")
            names.add(name)
            chrom = name.rsplit("_", 1)[0]
            start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
            attrs = dict(item.strip().split("=", 1) for item in parts[4].split(";") if "=" in item) if len(parts) > 4 else {}
            phase = int(attrs.get("phase", "0"))
            code = int(attrs.get("transl_table", attrs.get("genetic_code", "11")))
            partial = attrs.get("partial", "00")
            if start < 1 or end < start or strand not in (-1, 1) or phase not in (0, 1, 2):
                raise ValueError(f"Invalid coding coordinates for {name}")
            if code not in (1, 11) or len(partial) != 2 or set(partial) - {"0", "1"}:
                raise ValueError(f"Unsupported genetic code or partial flags for {name}")
            annotations.setdefault(chrom, []).append((name, start, end, strand, phase, code, partial))
    genes, reference = [], []
    seen = set()
    for header, sequence in fasta_records(reference_fasta):
        chrom = header.split()[0]
        if chrom in seen or chrom not in genomes:
            raise ValueError(f"Duplicate reference scaffold or missing STB entry: {chrom}")
        seen.add(chrom)
        sequence_id = hashlib.sha256(sequence.encode()).hexdigest()
        reference.append((chrom, genomes[chrom], len(sequence), sequence_id))
        for name, start, end, strand, phase, code, partial in annotations.pop(chrom, []):
            if end > len(sequence):
                raise ValueError(f"Gene {name} lies outside reference scaffold {chrom}")
            identity = [genomes[chrom], chrom, sequence_id, name, start, end, strand, phase, code]
            gene_id = int(content_digest(identity)[:16], 16)
            flags = partial if strand == 1 else partial[::-1]
            genes.append((gene_id, name, genomes[chrom], chrom, start, end, end-start+1,
                          strand, phase, code, flags[0] == "1", flags[1] == "1"))
    if annotations:
        raise ValueError("Annotated scaffold missing from reference: " + ", ".join(sorted(annotations)))
    if set(genomes) != seen:
        raise ValueError("STB contains scaffolds absent from reference")
    frame = pl.DataFrame(genes, schema=GENE_SCHEMA, orient="row").sort("gene_id")
    if frame["gene_id"].n_unique() != frame.height:
        raise ValueError("Gene ID collision")
    metadata = {REFERENCE_ID: content_digest(sorted(reference)),
                ANNOTATION_ID: content_digest(frame.rows()), SCHEMA_KEY: "2"}
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Stream reference codons in a second FASTA pass instead of retaining one
    # Python tuple for every coding triplet in a large reference collection.
    import pyarrow as pa
    import pyarrow.parquet as pq
    from zipstrain.utils import sha256_file
    codon_file = output.with_name("reference_codons.parquet")
    staging = codon_file.with_suffix(".part.parquet")
    schema = pa.schema([("gene_id", pa.uint64()), ("codon_index", pa.uint32()),
                        ("reference_codon", pa.string())]).with_metadata(
                            {k.encode():v.encode() for k,v in metadata.items()})
    by_chrom = frame.partition_by("chrom", as_dict=True)
    with pq.ParquetWriter(staging, schema, compression="zstd") as writer:
        batch = []
        for header, sequence in fasta_records(reference_fasta):
            chrom = header.split()[0]
            for gene in by_chrom.get((chrom,), pl.DataFrame(schema=GENE_SCHEMA)).iter_rows(named=True):
                coding = sequence[gene["start"]-1:gene["end"]]
                if gene["strand"] == -1:
                    coding = coding.translate(str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN"))[::-1]
                coding = coding[gene["phase"]:]
                for index in range(0, len(coding)-2, 3):
                    batch.append({"gene_id":gene["gene_id"], "codon_index":index//3,
                                  "reference_codon":coding[index:index+3]})
                    if len(batch) >= 65536:
                        writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                        batch.clear()
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
    staging.replace(codon_file)
    metadata["zipstrain_reference_codons_sha256"] = sha256_file(codon_file)
    frame.write_parquet(output, compression="zstd", metadata=metadata)
    return output


def read_gene_info(path: str | Path) -> pl.DataFrame:
    """Read and authenticate a V2 catalog; reject edited mappings and duplicate IDs."""
    path = Path(path)
    try:
        frame = pl.read_parquet(path).select(list(GENE_SCHEMA)).cast(GENE_SCHEMA).sort("gene_id")
    except pl.exceptions.PolarsError as exc:
        raise ValueError("V2 annotations require gene_info.parquet from prepare_profiling; rebuild legacy TSV annotations") from exc
    metadata = pl.read_parquet_metadata(path)
    if metadata.get(SCHEMA_KEY) != "2":
        raise ValueError("Unsupported gene catalog schema version")
    if not metadata.get(REFERENCE_ID) or metadata.get(ANNOTATION_ID) != content_digest(frame.rows()):
        raise ValueError("Gene annotation identity is missing or does not match its contents")
    if frame["gene_id"].n_unique() != frame.height:
        raise ValueError("Duplicate gene IDs")
    invalid = frame.filter((pl.col("start") < 1) | (pl.col("end") < pl.col("start")) |
        (pl.col("length") != pl.col("end")-pl.col("start")+1) |
        ~pl.col("strand").is_in([-1,1]) | ~pl.col("phase").is_in([0,1,2]) |
        ~pl.col("genetic_code").is_in([1,11]) |
        pl.any_horizontal([pl.col(name).is_null() for name in GENE_SCHEMA]))
    if invalid.height or frame["gene"].n_unique() != frame.height:
        raise ValueError("Invalid or unsupported gene annotation")
    return frame


def validate_profile_assets(profiles, gene_info: str | Path) -> dict[str, str]:
    """Validate sample/reference/annotation identities and return reusable metadata.

    Paths may point to position profiles or observed-codon tables. Missing
    identities are rejected rather than guessed from filenames or integer IDs.
    """
    read_gene_info(gene_info)
    expected = pl.read_parquet_metadata(gene_info)
    for path in profiles:
        actual = pl.read_parquet_metadata(path)
        if actual.get(SCHEMA_KEY) != "2":
            raise ValueError(f"Unsupported or missing V2 schema in {path}")
        for key in (REFERENCE_ID, ANNOTATION_ID):
            if actual.get(key) != expected[key]:
                raise ValueError(f"Incompatible {key} in {path}")
    return {key: expected[key] for key in (REFERENCE_ID, ANNOTATION_ID, SCHEMA_KEY)}


def read_gene_ranges(path) -> pl.LazyFrame:
    """Read new catalogs or legacy TSV intervals as ``gene, scaffold, start, end``.

    This adapter is for legacy matrix readers only; V2 comparisons authenticate
    the complete catalog with :func:`read_gene_info` instead.
    """
    path = Path(path).resolve()
    if path.suffix.lower() == ".parquet":
        return read_gene_info(path).lazy().select("gene", pl.col("chrom").alias("scaffold"), "start", "end")
    return pl.scan_csv(path, has_header=False, separator="\t",
                       new_columns=["gene", "scaffold", "start", "end"])


def range_sums(loci: pl.LazyFrame, genes: pl.DataFrame, columns, *, coordinates_sorted=False) -> pl.LazyFrame:
    """Sum numeric position columns within every gene, including overlaps.

    Input positions must be unique within each scaffold. Prefix sums are joined
    at inclusive ``end`` and ``start-1`` boundaries, avoiding per-position gene
    expansion. Uncovered genes receive zero; scaffold prefixes never mix.
    Set coordinates_sorted only for loci already ordered by chrom and pos.
    """
    columns = list(columns)
    ordered = loci if coordinates_sorted else loci.sort(["chrom", "pos"])
    prefix = ordered.with_columns(
        [pl.col(name).cast(pl.Float64).cum_sum().over("chrom").alias(name) for name in columns])
    bounds = pl.concat([
        genes.select("gene_id", "chrom", pl.col("end").alias("boundary")).with_columns(pl.lit(1).alias("sign")),
        genes.select("gene_id", "chrom", (pl.col("start")-1).alias("boundary")).with_columns(pl.lit(-1).alias("sign")),
    ]).lazy().sort(["chrom", "boundary"])
    return bounds.join_asof(prefix.select("chrom", "pos", *columns),
        left_on="boundary", right_on="pos", by="chrom", strategy="backward",
        check_sortedness=False).group_by("gene_id").agg(
            [(pl.col(name).fill_null(0)*pl.col("sign")).sum().alias(name) for name in columns])


def gene_statistics(profile: pl.LazyFrame, gene_info, *, min_cov=5, selection=None) -> pl.LazyFrame:
    """Build the agreed 15-column sample gene table from intervals and counts.

    Coverage and breadth use retained depth > 0, reference ANI uses known
    reference bases with depth >= ``min_cov``. Missing reference evidence and
    unavailable selection metrics are null, including undefined ratios.
    """
    if min_cov < 1:
        raise ValueError("min_cov must be positive")
    genes = read_gene_info(gene_info)
    depth = pl.sum_horizontal([pl.col(b).cast(pl.UInt64) for b in "ACGT"])
    has_reference = "ref_base_bitmask" in profile.collect_schema().names()
    ref = pl.col("ref_base_bitmask") if has_reference else pl.lit(0)
    callable_site = (depth >= min_cov) & ref.is_in([1, 2, 4, 8])
    match = pl.any_horizontal([(ref == bit) & (pl.col(base) > 0) for base, bit in zip("ACGT", [1, 2, 4, 8])])
    loci = profile.select("chrom", "pos", depth.alias("depth"), (depth > 0).alias("covered"),
                          callable_site.alias("callable"), (callable_site & match).alias("matches"))
    result = genes.lazy().join(range_sums(loci, genes, ["depth", "covered", "callable", "matches"]), on="gene_id", how="left").with_columns(
        (pl.col("depth").fill_null(0)/pl.col("length")).alias("coverage"),
        (pl.col("covered").fill_null(0)/pl.col("length")).alias("breadth"),
        pl.when(pl.col("callable") > 0).then(100*pl.col("matches")/pl.col("callable")).otherwise(None).alias("ref_ani"))
    if selection is not None:
        result = result.join(selection.select("gene_id", *SELECTION_COLUMNS), on="gene_id", how="left")
    else:
        result = result.with_columns([pl.lit(None, dtype=pl.Float64).alias(name) for name in SELECTION_COLUMNS])
    return result.select(GENE_STATS_COLUMNS)
