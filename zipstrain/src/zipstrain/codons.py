"""Read-observed codon evidence and explicit uncorrected selection summaries.

This module never reconstructs haplotypes by multiplying marginal frequencies.
See ``docs/v2.md`` for filters, estimators, exclusions, and runnable examples.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import heapq
import json
import re
import subprocess

import numpy as np
import polars as pl

from zipstrain import gene_info as gi, options


# NCBI tables 1 and 11 have identical internal amino-acid assignments.
_BASE_ORDER = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = dict(zip(("".join(c) for c in product(_BASE_ORDER, repeat=3)), _AMINO_ACIDS))
COUNT_SCHEMA = {"gene_id": pl.UInt64, "codon_index": pl.UInt32,
                "codon": pl.String, "count": pl.UInt64}
METHOD = "sns_snv_reference_opportunities_v2"
ALLELE_INTEGRATION_MODES = ("consensus", "weighted")
# Matched codons reduced at a time. The worst case of 4,096 observed codon pairs per locus
# needs about 1 MiB per row, so this bounds a pathological batch near 8 GiB and a typical one
# (a few pairs per locus) near 10 MiB. Results do not depend on it.
SELECTION_BATCH_ROWS = 8192


@lru_cache(maxsize=64)
def synonymous_sites(codon: str) -> float:
    """Return NG-style synonymous opportunities (0..3) for a sense codon.

    Each of the nine possible single-base mutations has weight 1/3. Invalid
    codons and stops are rejected. Tables 1 and 11 share this assignment.
    """
    if codon not in GENETIC_CODE or GENETIC_CODE[codon] == "*":
        raise ValueError(f"Expected a sense codon, received {codon!r}")
    return sum(GENETIC_CODE[codon[:i]+b+codon[i+1:]] == GENETIC_CODE[codon]
               for i in range(3) for b in "ACGT" if b != codon[i])/3.0


@lru_cache(maxsize=4096)
def codon_changes(source: str, target: str) -> tuple[float, float]:
    """Count synonymous/nonsynonymous changes across shortest sense-only paths.

    Average every valid nucleotide-change order equally, using the full source
    and target codons, including shared non-reference bases. Reject invalid or
    stop endpoints and pairs with no sense-only path. Identical codons return
    zero changes. This is a counting component, not a corrected rate estimator.
    """
    synonymous_sites(source)
    synonymous_sites(target)
    differences = [i for i in range(3) if source[i] != target[i]]
    paths = []
    for order in permutations(differences):
        current, synonymous = source, 0
        for index in order:
            mutated = current[:index]+target[index]+current[index+1:]
            if GENETIC_CODE[mutated] == "*":
                break
            synonymous += GENETIC_CODE[current] == GENETIC_CODE[mutated]
            current = mutated
        else:
            paths.append(synonymous)
    if not paths:
        raise ValueError("No sense-only substitution path")
    s = sum(paths)/len(paths)
    return s, len(differences)-s


def aligned_read_bases(position: int, cigar: str, sequence: str, qualities: str) -> dict:
    """Map one SAM alignment to 1-based reference bases and query indices.

    Return ``{position: (base, phred_quality, query_index)}``. Clipping,
    insertions, deletions, skips, padding, '=' and 'X' are honored. Reject
    malformed CIGARs, mismatched lengths, and missing qualities. Query indices
    let callers exclude triplets interrupted by an insertion.
    """
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if position < 1:
        raise ValueError("SAM reference positions must be one-based")
    if not tokens or "".join(n+op for n, op in tokens) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    if qualities == "*" or len(qualities) != len(sequence):
        raise ValueError("Missing or mismatched base qualities")
    query, reference, result = 0, position, {}
    for length_text, op in tokens:
        length = int(length_text)
        if length < 1:
            raise ValueError("CIGAR operations must be positive")
        if op in "M=X":
            if query+length > len(sequence):
                raise ValueError("CIGAR exceeds sequence length")
            for offset in range(length):
                result[reference+offset] = (sequence[query+offset].upper(), ord(qualities[query+offset])-33, query+offset)
            reference += length
            query += length
        elif op in "IS":
            query += length
        elif op in "DN":
            reference += length
    if query != len(sequence):
        raise ValueError("CIGAR query length does not match sequence")
    return result


def read_codon_observations(bases: dict, gene: dict, min_baseq: int = 13):
    """Yield ``(codon_index, coding_triplet, minimum_quality)`` for a read/gene.

    Require three consecutive query bases covering a complete reference codon.
    Exclude ambiguous bases, low-quality calls, insertions within a triplet,
    deletions, and skips. Minus-strand genes are reverse-complemented. Partial
    terminal codons and the explicit initial phase are excluded.
    """
    if not bases:
        return
    origin = gene["start"]+gene["phase"] if gene["strand"] == 1 else gene["end"]-gene["phase"]
    lo, hi = min(bases), max(bases)
    for index in range(max(0, (hi-origin)//-3) if gene["strand"] == -1 else max(0, (lo-origin)//3),
                       (gene["length"]-gene["phase"])//3):
        positions = [origin + gene["strand"]*(index*3+offset) for offset in range(3)]
        if gene["strand"] == 1 and positions[0] > hi:
            break
        if gene["strand"] == -1 and positions[0] < lo:
            break
        if not all(p in bases for p in positions):
            continue
        calls = [bases[p] for p in positions]
        if any(b not in "ACGT" or q < min_baseq for b, q, _ in calls):
            continue
        if any(calls[i+1][2]-calls[i][2] != gene["strand"] for i in (0, 1)):
            continue
        codon = "".join(b for b, _, _ in calls)
        if gene["strand"] == -1:
            codon = codon.translate(str.maketrans("ACGT", "TGCA"))
        yield index, codon, min(q for _, q, _ in calls)


def _extract_raw_codons(bam_file, genes: pl.DataFrame, raw, *, min_baseq, min_mapq, min_read_ani,
                        read_inclusion) -> None:
    """Write per-fragment codon counts for a coordinate-sorted BAM, read through samtools."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from zipstrain.profile import _sam_alignment_passes_profile_filters

    bins = {}
    for gene in genes.iter_rows(named=True):
        for bucket in range(gene["start"]//16384, gene["end"]//16384+1):
            bins.setdefault((gene["chrom"], bucket), []).append(gene)
    schema = pa.schema([("gene_id", pa.uint64()), ("codon_index", pa.uint32()),
                        ("codon", pa.string()), ("count", pa.uint64())])
    counts, active, expiry = Counter(), {}, []
    chrom, previous, finished = None, -1, set()
    process = subprocess.Popen(["samtools", "view", str(bam_file)], stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, text=True)
    try:
        with pq.ParquetWriter(raw, schema, compression="zstd") as writer:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                flag, start, mapq = int(fields[1]), int(fields[3]), int(fields[4])
                if flag & 0xF04:
                    continue
                if not _sam_alignment_passes_profile_filters(flag=flag, mapq=mapq, cigar=fields[5],
                    optional_fields=fields[11:], min_mapq=min_mapq, min_read_ani=min_read_ani,
                    read_inclusion=read_inclusion, rname=fields[2], rnext=fields[6]):
                    continue
                if fields[2] != chrom:
                    if fields[2] in finished:
                        raise ValueError("Codon profiling requires coordinate-sorted BAM")
                    if chrom is not None:
                        finished.add(chrom)
                    if counts:
                        writer.write_table(pa.Table.from_pylist([dict(gene_id=g, codon_index=i, codon=c, count=n)
                            for (g, i, c), n in counts.items() if n > 0], schema=schema))
                    counts, active, expiry = Counter(), {}, []
                    chrom, previous = fields[2], -1
                if start < previous:
                    raise ValueError("Codon profiling requires coordinate-sorted BAM")
                previous = start
                while expiry and expiry[0][0] < start:
                    _end, key = heapq.heappop(expiry)
                    active.pop(key, None)
                if fields[10] == "*":
                    continue
                bases = aligned_read_bases(start, fields[5], fields[9], fields[10])
                if not bases:
                    continue
                candidates = {}
                for bucket in range(start//16384, max(bases)//16384+1):
                    for gene in bins.get((chrom, bucket), []):
                        candidates[gene["gene_id"]] = gene
                rg = next((f for f in fields[11:] if f.startswith("RG:Z:")), "")
                for gene in candidates.values():
                    for index, codon, quality in read_codon_observations(bases, gene, min_baseq):
                        key = (fields[0], rg, gene["gene_id"], index)
                        old = active.get(key)
                        if old is not None:
                            old_codon, old_quality = old
                            if quality < old_quality or (quality == old_quality and codon == old_codon):
                                continue
                            if old_codon is not None:
                                counts[gene["gene_id"], index, old_codon] -= 1
                            if quality == old_quality:
                                active[key] = (None, quality)
                                continue
                        else:
                            heapq.heappush(expiry, (max(bases), key))
                        counts[gene["gene_id"], index, codon] += 1
                        active[key] = (codon, quality)
            if counts:
                writer.write_table(pa.Table.from_pylist([dict(gene_id=g, codon_index=i, codon=c, count=n)
                    for (g, i, c), n in counts.items() if n > 0], schema=schema))
        if process.wait() != 0:
            raise RuntimeError("samtools view failed during codon profiling")
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        process.stdout.close()


def write_codon_counts(bam_file, profile_file, gene_info, output_file, *,
                       min_baseq=13, min_mapq=0, min_read_ani=options.PROFILE_MIN_READ_ANI,
                       read_inclusion=options.PROFILE_READ_INCLUSION, extractor=None) -> Path:
    """Collect actual full-read triplets from a coordinate-sorted BAM.

    Counts are per fragment: overlapping mates with the same QNAME/RG contribute
    once per gene/codon. Choose the higher minimum-quality triplet; equal-quality
    disagreement discards that observation. Never bridge an unobserved mate gap.
    The active mate cache expires as reference coordinates advance. Observations
    are flushed per scaffold, then gated against alleles retained in the position
    profile. Base qualities here are original SAM qualities, not BAQ-adjusted;
    the retained-position gate also applies the profiler's error/frequency policy.
    Secondary, supplementary, duplicate, QC-failed and unmapped reads are excluded.

    ``extractor`` replaces the Python extraction with a faster equivalent: a callable taking
    ``(bam_file, genes, raw_path)`` and the keyword filters. The gate and file are shared.
    """
    from zipstrain.utils import sha256_file

    if min_baseq < 0 or min_mapq < 0:
        raise ValueError("Quality thresholds must be nonnegative")
    metadata = gi.validate_profile_assets([profile_file], gene_info)
    genes = gi.read_gene_info(gene_info)
    metadata.update({"zipstrain_codon_evidence": "read_phased_v1",
                     "zipstrain_source_profile_sha256": sha256_file(profile_file),
                     "zipstrain_codon_quality": "original_phred_then_retained_position_gate",
                     "zipstrain_codon_filters": json.dumps({"min_baseq": min_baseq, "min_mapq": min_mapq,
                         "min_read_ani": min_read_ani, "read_inclusion": read_inclusion}, sort_keys=True),
                     "zipstrain_selection_method": METHOD})
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = output.with_suffix(".raw.parquet")
    try:
        (extractor or _extract_raw_codons)(bam_file, genes, raw, min_baseq=min_baseq, min_mapq=min_mapq,
                                           min_read_ani=min_read_ani, read_inclusion=read_inclusion)
        # A triplet must consist entirely of alleles retained by the profiler.
        observed = pl.scan_parquet(raw).cast(COUNT_SCHEMA).join(genes.lazy(), on="gene_id", how="left")
        valid = observed
        for offset in range(3):
            pos = pl.when(pl.col("strand") == 1).then(pl.col("start")+pl.col("phase")+3*pl.col("codon_index")+offset).otherwise(
                pl.col("end")-pl.col("phase")-3*pl.col("codon_index")-offset).cast(pl.Int64)
            base = pl.col("codon").str.slice(offset, 1)
            base = pl.when(pl.col("strand") == 1).then(base).otherwise(base.replace_strict(dict(zip("ACGT", "TGCA"))))
            valid = valid.with_columns(pos.alias("pos"), base.alias("_base")).join(
                pl.scan_parquet(profile_file).select("chrom", pl.col("pos").cast(pl.Int64), "A", "C", "G", "T"),
                on=["chrom", "pos"], how="inner").filter(
                    pl.any_horizontal([(pl.col("_base") == b) & (pl.col(b) > 0) for b in "ACGT"])).drop("A", "C", "G", "T", "_base", "pos")
        valid.select(list(COUNT_SCHEMA)).sort(["gene_id", "codon_index", "codon"]).sink_parquet(
            output, compression="zstd", metadata=metadata)
    finally:
        raw.unlink(missing_ok=True)
    return output


def read_codon_counts(path, gene_info) -> pl.LazyFrame:
    """Authenticate and validate observed-codon evidence before selection analysis.

    Reject marginal/reconstructed evidence, invalid triplets, zero counts,
    duplicate haplotype rows, unknown genes, and out-of-range codon indices.
    The adjacent reference-codon file must carry the same annotation identity.
    """
    gi.validate_profile_assets([path], gene_info)
    if pl.read_parquet_metadata(path).get("zipstrain_codon_evidence") != "read_phased_v1":
        raise ValueError("Selection analysis requires read-phased codon evidence")
    reference_file = Path(gene_info).resolve().with_name("reference_codons.parquet")
    gi.validate_profile_assets([reference_file], gene_info)
    from zipstrain.utils import sha256_file
    if sha256_file(reference_file) != pl.read_parquet_metadata(gene_info).get("zipstrain_reference_codons_sha256"):
        raise ValueError("Reference codon contents do not match the preprofile catalog")
    frame = pl.scan_parquet(path).select(list(COUNT_SCHEMA)).cast(COUNT_SCHEMA)
    if frame.filter(pl.col("count").is_null() | (pl.col("count") == 0) | pl.col("codon").is_null() |
                    ~pl.col("codon").str.contains(r"^[ACGT]{3}$")).limit(1).collect(engine="streaming").height:
        raise ValueError("Invalid observed codon or count")
    if _has_duplicate_rows(frame):
        raise ValueError("Duplicate observed codon rows")
    # The reference holds codons 0 .. (length - phase) // 3 - 1 of every catalog gene, and its
    # digest was checked above, so the small catalog stands in for the full reference table.
    codon_counts = gi.read_gene_info(gene_info).lazy().select(
        "gene_id", ((pl.col("length") - pl.col("phase")) // 3).alias("_codons"))
    if frame.join(codon_counts, on="gene_id", how="left").filter(
            pl.col("_codons").is_null() | (pl.col("codon_index") >= pl.col("_codons"))
    ).limit(1).collect(engine="streaming").height:
        raise ValueError("Unknown gene ID or codon index")
    return frame


def _has_duplicate_rows(frame: pl.LazyFrame) -> bool:
    """Return whether a gene, codon index and codon appear twice.

    ``write_codon_counts`` writes tables sorted by those columns, where a duplicate is a row
    equal to the one before it, checked over streamed batches. Other tables are hashed whole.
    """
    keys = ["gene_id", "codon_index", "codon"]
    same = [pl.col(key) == pl.col(key).shift(1) for key in keys]
    later = [pl.col(key) > pl.col(key).shift(1) for key in keys]
    ordered = later[0] | (same[0] & (later[1] | (same[1] & (later[2] | same[2]))))
    checks = [ordered.fill_null(True).all().alias("sorted"),
              (same[0] & same[1] & same[2]).fill_null(False).any().alias("duplicate")]
    previous = None
    for batch in frame.select(keys).collect_batches(engine="streaming"):
        if not batch.height:
            continue
        rows = batch if previous is None else pl.concat([previous, batch])
        is_sorted, duplicate = rows.select(checks).row(0)
        if not is_sorted:
            return frame.select(pl.struct(*keys).is_duplicated().any()).collect().item()
        if duplicate:
            return True
        previous = batch.tail(1)
    return False


@lru_cache(maxsize=12288)
def site_changes(source: str, target: str, offset: int) -> tuple[float, float]:
    """Average the effect of one differing site across full-codon sense-only paths.

    This preserves the SNS/SNV site definitions from ``gene-overlap`` while
    fixing its reference-background substitution bug. Other shared nonreference
    bases remain present; multiple differing sites use all shortest sense-only
    orders. The three per-site results sum to :func:`codon_changes`.
    """
    synonymous_sites(source)
    synonymous_sites(target)
    if offset not in (0, 1, 2):
        raise ValueError("Codon offset must be 0, 1, or 2")
    if source[offset] == target[offset]:
        return 0.0, 0.0
    differences = [i for i in range(3) if source[i] != target[i]]
    effects = []
    for order in permutations(differences):
        current, synonymous = source, None
        for index in order:
            mutated = current[:index]+target[index]+current[index+1:]
            if GENETIC_CODE[mutated] == "*":
                break
            if index == offset:
                synonymous = float(GENETIC_CODE[current] == GENETIC_CODE[mutated])
            current = mutated
        else:
            effects.append(synonymous)
    if not effects:
        raise ValueError("No sense-only substitution path")
    s = sum(effects)/len(effects)
    return s, 1.0-s


def sibling_codon_table(profile_path) -> Path | None:
    """Return where a profile's codon table sits by naming convention, if it has one."""
    path = Path(profile_path)
    if not path.name.endswith("_profile.parquet"):
        return None
    return path.with_name(path.name.replace("_profile.parquet", "_codon_counts.parquet"))


def selection_frame(totals: dict) -> pl.LazyFrame:
    """Normalize SNS/SNV change totals by callable reference opportunities.

    Each gene's totals are ``[SNS_N, SNS_S, SNV_N, SNV_S, S_sites, N_sites]``.
    Return only ``gene_id`` and the six agreed metrics. No multiple-hit correction
    is applied. Zero opportunity denominators and zero dS/pS ratios yield null.
    """
    gene_ids = np.fromiter(totals, dtype=np.uint64, count=len(totals))
    values = np.array(list(totals.values()), dtype=np.float64).reshape(-1, 6)
    return _selection_columns(gene_ids, values)


def _selection_columns(gene_ids: np.ndarray, totals: np.ndarray) -> pl.LazyFrame:
    """Return :func:`selection_frame` rows for genes and their ``(n, 6)`` totals array."""
    names = ["_sns_n", "_sns_s", "_snv_n", "_snv_s", "_s_sites", "_n_sites"]
    frame = pl.DataFrame({"gene_id": pl.Series(gene_ids, dtype=pl.UInt64),
                          **{name: totals[:, index] for index, name in enumerate(names)}})
    per_site = lambda count, sites: pl.when(pl.col(sites) != 0).then(pl.col(count) / pl.col(sites))
    ratio = lambda top, bottom: pl.when(pl.col(bottom) != 0).then(pl.col(top) / pl.col(bottom))
    return frame.lazy().select(
        "gene_id", per_site("_sns_n", "_n_sites").alias("dN"), per_site("_sns_s", "_s_sites").alias("dS"),
        per_site("_snv_n", "_n_sites").alias("pN"), per_site("_snv_s", "_s_sites").alias("pS"),
    ).with_columns(dN_dS=ratio("dN", "dS"), pN_pS=ratio("pN", "pS")).select("gene_id", *gi.SELECTION_COLUMNS)


_CODON_IDS = {codon: index for index, codon in enumerate(GENETIC_CODE)}
_BASES = np.array([["ACGT".index(base) for base in codon] for codon in _CODON_IDS], dtype=np.uint8)
_OPPORTUNITIES = {codon: synonymous_sites(codon) for codon, amino_acid in GENETIC_CODE.items() if amino_acid != "*"}
_KEYS = ["gene_id", "codon_index"]
# Codons packed or matched at a time, in ranges of whole genes. This bounds memory near a few
# hundred MB however large the catalog is. Results do not depend on it.
SELECTION_CHUNK_CODONS = 1_000_000
_NOT_SENSE = 255  # reference lookup entry for stop and ambiguous reference codons
_SITES_BY_ID = np.full(256, np.nan)
_SITES_BY_ID[[_CODON_IDS[codon] for codon in _OPPORTUNITIES]] = list(_OPPORTUNITIES.values())
_PACKED_SCHEMA = {"gene_id": pl.UInt64, "codon_index": pl.UInt32, "codon": pl.List(pl.UInt8),
                  "count": pl.List(pl.UInt64), "depth": pl.UInt64, "s_sites": pl.Float64}


def _selection_arguments(min_cov, min_freq, allele_integration, batch_size=1) -> None:
    """Reject selection settings that have no defined meaning."""
    if allele_integration not in ALLELE_INTEGRATION_MODES:
        raise ValueError("allele_integration must be consensus or weighted")
    if min_cov < 1 or not 0 <= min_freq <= 1:
        raise ValueError("Invalid codon coverage or frequency threshold")
    if batch_size < 1:
        raise ValueError("Selection batch size must be positive")


@lru_cache(maxsize=1)
def _consequence_lookup() -> np.ndarray:
    """Return nonsynonymous and synonymous effects by offset and ordered codon pair."""
    effects = np.zeros((3, 64, 64, 2), dtype=np.float64)
    for source in _OPPORTUNITIES:
        for target in _OPPORTUNITIES:
            for offset in range(3):
                if source[offset] != target[offset]:
                    try:
                        synonymous, nonsynonymous = site_changes(source, target, offset)
                    except ValueError:
                        continue
                    effects[offset, _CODON_IDS[source], _CODON_IDS[target]] = nonsynonymous, synonymous
    return effects


def _run_starts(values: np.ndarray) -> np.ndarray:
    """Return where each run of equal neighbouring values starts."""
    return np.flatnonzero(np.r_[True, values[1:] != values[:-1]]) if len(values) else np.zeros(0, dtype=np.intp)


def _run_positions(values: np.ndarray, table: np.ndarray) -> np.ndarray:
    """Return ``np.searchsorted(table, values)``, searching once per run of equal values."""
    starts = _run_starts(values)
    return np.repeat(np.searchsorted(table, values[starts]), np.diff(np.r_[starts, len(values)]))


class _ReferenceLookup:
    """Reference codon IDs of a set of genes, one byte per codon, by gene ID and codon index.

    ``reference_codons.parquet`` is in scaffold order, not gene-ID order, so it cannot be
    merged with sorted codon tables. One streaming read into a flat array replaces that join.
    """

    def __init__(self, catalog, genes: pl.DataFrame, *, subset: bool = False):
        """``genes`` is the catalog sorted by gene ID, or, with ``subset``, part of it."""
        self.subset = subset
        self.gene_ids = genes["gene_id"].to_numpy()
        sizes = ((genes["length"] - genes["phase"]) // 3).to_numpy().astype(np.int64)
        self.offsets = np.concatenate([[0], np.cumsum(sizes)])
        self.ids = np.full(self.offsets[-1], _NOT_SENSE, dtype=np.uint8)
        scan = pl.scan_parquet(Path(catalog).resolve().with_name("reference_codons.parquet"))
        if subset:
            scan = scan.join(genes.lazy().select("gene_id"), on="gene_id", how="semi")
        scan = scan.select(*_KEYS, pl.col("reference_codon").replace_strict(
            {codon: _CODON_IDS[codon] for codon in _OPPORTUNITIES},
            default=_NOT_SENSE, return_dtype=pl.UInt8).alias("id"))
        for batch in scan.collect_batches(engine="streaming"):
            self.ids[self._index(batch)] = batch["id"].to_numpy()

    @classmethod
    def for_genome(cls, catalog, genome: str = "all") -> "_ReferenceLookup":
        """Return the lookup of every catalog gene, or of one genome's genes."""
        genes = gi.read_gene_info(catalog)
        if genome == "all":
            return cls(catalog, genes)
        return cls(catalog, genes.filter(pl.col("genome") == genome), subset=True)

    def _index(self, frame: pl.DataFrame) -> np.ndarray:
        where = _run_positions(frame["gene_id"].to_numpy(), self.gene_ids)
        return self.offsets[where] + frame["codon_index"].to_numpy()

    def codon_ids(self, frame: pl.DataFrame) -> np.ndarray:
        """Return the reference codon ID of each row, or ``_NOT_SENSE``."""
        return self.ids[self._index(frame)]

    def gene_ranges(self, chunk: int) -> list[tuple[int, int]]:
        """Split the genes into consecutive (first, last) gene-ID ranges of about ``chunk`` codons."""
        if not len(self.gene_ids):
            return []
        group = self.offsets[:-1] // chunk
        starts = np.flatnonzero(np.diff(group, prepend=-1))
        ends = np.append(starts[1:], len(self.gene_ids)) - 1
        return [(int(self.gene_ids[s]), int(self.gene_ids[e])) for s, e in zip(starts, ends)]


def _gene_range(first: int | None, last: int | None) -> pl.Expr:
    """Rows with ``first < gene_id <= last``; either bound may be open."""
    where = pl.lit(True)
    if first is not None:
        where &= pl.col("gene_id") > pl.lit(first, dtype=pl.UInt64)
    if last is not None:
        where &= pl.col("gene_id") <= pl.lit(last, dtype=pl.UInt64)
    return where


def _pack_callable_codons(source, catalog, destination, min_cov: int, genome: str = "all",
                          metadata: dict[str, str] | None = None, min_freq: float = 0.0) -> None:
    """Stage one sample as sorted sparse codon groups after validating its evidence.

    Each row is one callable codon: its observed codon IDs and counts, retained depth, and
    the reference codon's synonymous opportunities. ``min_freq`` is an optional additional
    *codon* frequency floor; position-level allele filtering has already occurred, so it is
    zero by default. Reference stop or ambiguous codons are excluded. Genes are packed in
    gene-ID ranges of about ``SELECTION_CHUNK_CODONS`` codons, so memory does not grow with the table.
    """
    import pyarrow.parquet as pq

    counts = read_codon_counts(source, catalog)
    reference = _ReferenceLookup.for_genome(catalog, genome)
    schema = pl.DataFrame(schema=_PACKED_SCHEMA).to_arrow().schema.with_metadata(metadata or {})
    # Packed files are short-lived staging, so favour write speed over size.
    with pq.ParquetWriter(destination, schema, compression="zstd", compression_level=1) as writer:
        for part in _packed_ranges(counts, reference, min_cov, min_freq):
            writer.write_table(part.to_arrow().cast(schema))


def _packed_ranges(counts: pl.LazyFrame, reference: _ReferenceLookup, min_cov: int, min_freq: float = 0.0):
    """Yield one sample's packed callable codons, a range of whole genes at a time.

    ``counts`` is a validated codon table and ``reference`` the lookup of the genes to pack.
    Two samples packed with the same lookup yield the same ranges, in gene order.
    """
    for first, last in reference.gene_ranges(SELECTION_CHUNK_CODONS):
        part = counts.filter(pl.col("gene_id").is_between(pl.lit(first, dtype=pl.UInt64),
                                                          pl.lit(last, dtype=pl.UInt64)))
        if reference.subset:
            part = part.filter(pl.col("gene_id").is_in(pl.Series(reference.gene_ids).implode()))
        if min_freq > 0:
            reads = pl.col("count").cast(pl.Float64)
            part = part.filter(reads >= reads.sum().over(_KEYS) * min_freq)
        part = _group_codons(part.with_columns(
            pl.col("codon").replace_strict(_CODON_IDS, return_dtype=pl.UInt8)).collect())
        part = part.filter(pl.col("depth") >= min_cov)
        sites = _SITES_BY_ID[reference.codon_ids(part)]
        yield part.with_columns(pl.Series("s_sites", sites, nan_to_null=True)).drop_nulls("s_sites") \
            .cast(_PACKED_SCHEMA)


def _group_codons(rows: pl.DataFrame) -> pl.DataFrame:
    """Group observed codon rows into one row per locus, sorted by gene and codon index.

    Return each locus's codon IDs and counts, in input order, and their summed depth.
    Tables from ``write_codon_counts`` are sorted, so loci are runs of rows, found and summed
    in NumPy and wrapped as list columns without copying. Other rows are sorted first.
    """
    import pyarrow as pa

    gene, index = rows["gene_id"].to_numpy(), rows["codon_index"].to_numpy()
    if not np.all((gene[1:] > gene[:-1]) | ((gene[1:] == gene[:-1]) & (index[1:] >= index[:-1]))):
        rows = rows.sort(_KEYS, maintain_order=True)
        gene, index = rows["gene_id"].to_numpy(), rows["codon_index"].to_numpy()
    starts = np.flatnonzero(np.r_[True, (gene[1:] != gene[:-1]) | (index[1:] != index[:-1])]) \
        if len(gene) else np.zeros(0, dtype=np.int64)
    reads = rows["count"].to_numpy()
    lengths = np.diff(np.r_[starts, len(gene)])
    if len(reads) and int(reads.max()) * int(lengths.max()) >= 2**64:  # sum exactly; a too-large depth fails the cast
        return (rows.group_by(_KEYS, maintain_order=True)
                .agg(pl.col("codon"), pl.col("count"), pl.col("count").cast(pl.Int128).sum().alias("depth"))
                .with_columns(pl.col("depth").cast(pl.UInt64)))
    offsets = pa.array(np.r_[starts, len(gene)].astype(np.int64))
    as_list = lambda name: pl.from_arrow(pa.LargeListArray.from_arrays(offsets, rows[name].to_arrow()))
    return pl.DataFrame({
        "gene_id": gene[starts], "codon_index": index[starts],
        "codon": as_list("codon"), "count": as_list("count"),
        "depth": np.add.reduceat(reads, starts) if len(starts) else np.zeros(0, dtype=np.uint64),
    }, schema={key: _PACKED_SCHEMA[key] for key in ("gene_id", "codon_index", "codon", "count", "depth")})


def _file_gene_ranges(path, chunk_codons: int) -> list[pl.Expr]:
    """Split a gene-sorted parquet file into filters of about ``chunk_codons`` rows of whole genes.

    Cuts fall at the last gene ID of a row group, so a gene split across row groups stays in
    one range, and each filter reads only the row groups it needs. A file without row-group
    statistics is one range.
    """
    import pyarrow.parquet as pq

    metadata = pq.ParquetFile(path).metadata
    column = metadata.schema.names.index("gene_id")
    cuts, rows = [], 0
    for index in range(metadata.num_row_groups):
        group = metadata.row_group(index)
        statistics = group.column(column).statistics
        if statistics is None or not statistics.has_min_max:
            cuts = []
            break
        rows += group.num_rows
        if rows >= chunk_codons:
            cuts.append(statistics.max)
            rows = 0
    return [_gene_range(first, last) for first, last in zip([None, *cuts], [*cuts, None])]


def _batch_totals(batch: pl.DataFrame, integration: str) -> np.ndarray:
    """Reduce a block of matched codons to six numerator and opportunity totals per codon."""
    left = batch["codon"].explode().to_numpy()
    right = batch["codon_2"].explode().to_numpy()
    reads_a = batch["count"].explode().to_numpy()
    reads_b = batch["count_2"].explode().to_numpy()
    n_a = batch["codon"].list.len().cast(pl.Int64).to_numpy()
    n_b = batch["codon_2"].list.len().cast(pl.Int64).to_numpy()
    depth_a = batch["depth"].to_numpy()
    depth_b = batch["depth_2"].to_numpy()
    n = batch.height
    rows = np.arange(n, dtype=np.intp)
    row_a, row_b = np.repeat(rows, n_a), np.repeat(rows, n_b)
    combinations_per_row = n_a * n_b
    row_pair = np.repeat(rows, combinations_per_row)
    local = np.arange(len(row_pair)) - np.repeat(
        np.cumsum(combinations_per_row) - combinations_per_row, combinations_per_row
    )
    index_a = (np.cumsum(n_a) - n_a)[row_pair] + local // n_b[row_pair]
    index_b = (np.cumsum(n_b) - n_b)[row_pair] + local % n_b[row_pair]
    codon_a, codon_b = left[index_a], right[index_b]

    counts_product = reads_a[index_a].astype(np.float64) * reads_b[index_b]
    if integration == "weighted":
        weights = counts_product / (depth_a[row_pair].astype(np.float64) * depth_b[row_pair])

    totals = np.zeros((n, 6), dtype=np.float64)
    totals[:, 4] = batch["s_sites"].to_numpy()
    totals[:, 5] = 3 - totals[:, 4]
    effects = _consequence_lookup()
    for offset in range(3):
        a = np.zeros((n, 4), dtype=np.uint64)
        b = np.zeros_like(a)
        np.add.at(a, (row_a, _BASES[left, offset]), reads_a)
        np.add.at(b, (row_b, _BASES[right, offset]), reads_b)
        fixed = ((a > 0).sum(axis=1) == 1) & ((b > 0).sum(axis=1) == 1)
        if integration == "consensus":
            # Support is the exact product of two integer allele counts. Use
            # Python integers only for unusual inputs that can overflow UInt64.
            if int(a.max()) * int(b.max()) > np.iinfo(np.uint64).max:
                support = a.astype(object)[:, :, None] * b.astype(object)[:, None, :]
            else:
                support = a[:, :, None] * b[:, None, :]
            support[:, np.arange(4), np.arange(4)] = 0
            top = support.max(axis=(1, 2), keepdims=True)
            best = (support == top) & (support > 0)
            ties = best.sum(axis=(1, 2))
            index = (row_pair, _BASES[codon_a, offset], _BASES[codon_b, offset])
            chosen = best[index]
            denominator = support[index].astype(np.float64) * ties[row_pair]
            weights = np.divide(counts_product, denominator, out=np.zeros_like(counts_product), where=chosen)
        for effect in range(2):
            sums = np.bincount(row_pair, weights=weights * effects[offset, codon_a, codon_b, effect], minlength=n)
            totals[:, effect] += np.where(fixed, sums, 0)
            totals[:, effect + 2] += np.where(fixed, 0, sums)
    return totals


def _keys_sorted(frame: pl.DataFrame) -> bool:
    """Return whether rows are in strictly increasing (gene_id, codon_index) order."""
    gene, index = frame["gene_id"].to_numpy(), frame["codon_index"].to_numpy()
    return bool(np.all((gene[1:] > gene[:-1]) | ((gene[1:] == gene[:-1]) & (index[1:] > index[:-1]))))


def _single_codon(frame: pl.DataFrame) -> np.ndarray:
    """Return each row's codon ID when it holds exactly one observed codon, else -1."""
    codon = pl.col("codon")
    return frame.select(pl.when(codon.list.len() == 1).then(codon.list.first().cast(pl.Int16))
                        .otherwise(-1))["codon"].to_numpy()


def _match_sorted(a: pl.DataFrame, b: pl.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return the rows of ``a`` and ``b`` holding the same codon, in ``a``'s order.

    Both frames are sorted by gene and codon index. Each locus becomes one integer key, the
    gene's rank in this range and the codon index, so Polars can merge two sorted columns.
    """
    gene_a, gene_b = a["gene_id"].to_numpy(), b["gene_id"].to_numpy()
    genes = np.union1d(gene_a[_run_starts(gene_a)], gene_b[_run_starts(gene_b)])

    def keys(gene, frame, name):
        key = (_run_positions(gene, genes).astype(np.uint64) << np.uint64(32)) | \
            frame["codon_index"].to_numpy().astype(np.uint64)
        return pl.DataFrame({"key": key}).with_row_index(name).set_sorted("key")

    pairs = keys(gene_a, a, "i").join(keys(gene_b, b, "j"), on="key", maintain_order="left")
    return pairs["i"].to_numpy(), pairs["j"].to_numpy()


def _selection_totals(matched, integration: str, batch_size: int) -> pl.LazyFrame:
    """Reduce matched codons, range by range, to one row of selection metrics per gene.

    ``matched`` yields ``(a, b, i, j)``: two samples' codons for a range of whole genes, sorted
    by gene and codon index, and the rows that match, in ``a``'s order. ``a`` carries the
    reference opportunities. Opportunities are integer multiples of 1/3, so they are summed
    exactly. A codon that is one and the same observed codon in both samples changes nothing,
    so only the other rows are reduced, ``batch_size`` at a time; each is added to its gene's
    totals in row order, so the sums depend on neither ``batch_size`` nor the ranges.
    """
    gene_ids, totals = [np.zeros(0, dtype=np.uint64)], [np.zeros((0, 6))]
    for a, b, i, j in matched:
        if not len(i):
            continue
        genes = a["gene_id"].to_numpy()[i]
        starts = _run_starts(genes)
        rank = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, len(genes)]))
        range_totals = np.zeros((len(starts), 6))
        thirds = np.rint(a["s_sites"].to_numpy()[i] * 3)
        synonymous = np.bincount(rank, weights=thirds, minlength=len(starts))
        loci = np.bincount(rank, minlength=len(starts))
        range_totals[:, 4] = synonymous / 3
        range_totals[:, 5] = (9 * loci - synonymous) / 3
        single_a, single_b = _single_codon(a)[i], _single_codon(b)[j]
        changing = np.flatnonzero((single_a != single_b) | (single_a < 0))
        columns = ["codon", "count", "depth"]
        for start in range(0, len(changing), batch_size):
            rows = changing[start:start + batch_size]
            batch = a[i[rows]].select("gene_id", "s_sites", *columns).hstack(
                b[j[rows]].select(pl.col(columns).name.suffix("_2")))
            np.add.at(range_totals[:, :4], rank[rows], _batch_totals(batch, integration)[:, :4])
        gene_ids.append(genes[starts])
        totals.append(range_totals)
    return _selection_columns(np.concatenate(gene_ids), np.concatenate(totals))


def vectorized_pair_selection(left, right, *, allele_integration: str = options.ALLELE_INTEGRATION,
                              batch_size: int = SELECTION_BATCH_ROWS) -> pl.LazyFrame:
    """Compare two packed samples one gene range at a time, in bounded batches.

    Results depend on neither ``batch_size`` nor ``SELECTION_CHUNK_CODONS``.
    """
    _selection_arguments(1, 0.0, allele_integration, batch_size)
    first, second = pl.scan_parquet(left), pl.scan_parquet(right).drop("s_sites")

    def matched():
        for where in _file_gene_ranges(left, SELECTION_CHUNK_CODONS):
            a, b = pl.collect_all([first.filter(where), second.filter(where)])
            a, b = (frame if _keys_sorted(frame) else frame.sort(_KEYS) for frame in (a, b))
            yield a, b, *_match_sorted(a, b)

    return _selection_totals(matched(), allele_integration, batch_size)


def reference_selection(codon_counts, gene_info, *, min_cov=options.COMPARE_MIN_COV, min_freq=0.0,
                        allele_integration=options.ALLELE_INTEGRATION,
                        batch_size=SELECTION_BATCH_ROWS) -> pl.LazyFrame:
    """Compute reference-relative SNS dN/dS and SNV pN/pS using the branch method.

    Callability and opportunities are per complete reference codon. Classification
    is per nucleotide site, so one polymorphic codon can contribute both fixed
    SNS and SNV events at different offsets. Consensus selects the most-supported
    nonreference allele (exact ties averaged); weighted mode uses all alternatives
    at their observed frequency. Consequences use full reference/observed codon
    paths conditional on the chosen allele, with contexts weighted by support.
    Invalid/stop consequences are excluded from change counts. Reference stop or
    ambiguous codons are excluded from opportunities. No Jukes–Cantor correction.

    This is the pair estimator with the reference as a sample holding one read of
    its own codon at every locus.
    """
    _selection_arguments(min_cov, min_freq, allele_integration, batch_size)
    counts = read_codon_counts(codon_counts, gene_info)
    reference = _ReferenceLookup.for_genome(gene_info)

    def matched():
        # Every packed codon has a sense reference codon, which sample 1 holds one read of.
        for b in _packed_ranges(counts, reference, min_cov, min_freq):
            a = b.select(*_KEYS, "s_sites").with_columns(
                pl.Series("reference_id", reference.codon_ids(b), dtype=pl.UInt8)
            ).select(*_KEYS, "s_sites", pl.concat_list("reference_id").alias("codon"),
                     pl.concat_list(pl.lit(1, dtype=pl.UInt64)).alias("count"),
                     pl.lit(1, dtype=pl.UInt64).alias("depth"))
            rows = np.arange(b.height)
            yield a, b, rows, rows

    return _selection_totals(matched(), allele_integration, batch_size)


def pairwise_selection(codon_counts_1, codon_counts_2, gene_info, *, min_cov=options.COMPARE_MIN_COV,
                       min_freq=0.0, allele_integration=options.ALLELE_INTEGRATION,
                       batch_size=SELECTION_BATCH_ROWS) -> pl.LazyFrame:
    """Compute pairwise SNS dN/dS and SNV pN/pS with the branch's integration rules.

    Merge sorted callable codons from both samples. Use reference opportunities
    over their intersection. Consensus selects the most-supported differing
    allele pair at each polymorphic site, averaging exact ties; weighted mode
    uses all differing pairs at their independent observed frequencies. Phased
    contexts are integrated conditionally, with full-codon path averaging. Thus
    GGA/GGT stays synonymous even against reference AAA. Identical polymorphic
    populations can have positive pN/pS: this measures polymorphic allele-pair
    consequences, not only substitutions separating populations. Raw proportions,
    no Jukes–Cantor correction; six output metrics, undefined ratios null.
    """
    _selection_arguments(min_cov, min_freq, allele_integration, batch_size)
    first, second = read_codon_counts(codon_counts_1, gene_info), read_codon_counts(codon_counts_2, gene_info)
    reference = _ReferenceLookup.for_genome(gene_info)

    def matched():
        # Both samples are packed in the same gene ranges, so range k of one meets range k of the other.
        for a, b in zip(_packed_ranges(first, reference, min_cov, min_freq),
                        _packed_ranges(second, reference, min_cov, min_freq), strict=True):
            b = b.drop("s_sites")
            yield a, b, *_match_sorted(a, b)

    return _selection_totals(matched(), allele_integration, batch_size)
