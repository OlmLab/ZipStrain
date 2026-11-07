"""Read-observed codon evidence and explicit uncorrected selection summaries.

This module never reconstructs haplotypes by multiplying marginal frequencies.
See ``docs/v2.md`` for filters, estimators, exclusions, and runnable examples.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
from tempfile import TemporaryDirectory
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
                    ~pl.col("codon").str.contains(r"^[ACGT]{3}$")).limit(1).collect().height:
        raise ValueError("Invalid observed codon or count")
    if frame.select(pl.struct("gene_id", "codon_index", "codon").is_duplicated().any()).collect().item():
        raise ValueError("Duplicate observed codon rows")
    reference = pl.scan_parquet(reference_file).select("gene_id", "codon_index")
    if frame.join(reference, on=["gene_id", "codon_index"], how="anti").limit(1).collect().height:
        raise ValueError("Unknown gene ID or codon index")
    return frame


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
    rows = []
    for gene_id, (dn, ds, pn, ps, s_sites, n_sites) in totals.items():
        d_n, d_s = dn/n_sites if n_sites else None, ds/s_sites if s_sites else None
        p_n, p_s = pn/n_sites if n_sites else None, ps/s_sites if s_sites else None
        rows.append((gene_id, d_n, d_s, d_n/d_s if d_s and d_n is not None else None,
                     p_n, p_s, p_n/p_s if p_s and p_n is not None else None))
    return pl.DataFrame(rows, schema={"gene_id": pl.UInt64, **{c: pl.Float64 for c in gi.SELECTION_COLUMNS}}, orient="row").lazy()


_CODON_IDS = {codon: index for index, codon in enumerate(GENETIC_CODE)}
_BASES = np.array([["ACGT".index(base) for base in codon] for codon in _CODON_IDS], dtype=np.uint8)
_OPPORTUNITIES = {codon: synonymous_sites(codon) for codon, amino_acid in GENETIC_CODE.items() if amino_acid != "*"}
_KEYS = ["gene_id", "codon_index"]
_NO_TOTALS = (0.0,) * 6


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


def _pack_callable_codons(source, catalog, destination, min_cov: int, genome: str = "all",
                          metadata: dict[str, str] | None = None, min_freq: float = 0.0) -> None:
    """Stage one sample as sorted sparse codon groups after validating its evidence.

    Each row is one callable codon: its observed codon IDs and counts, retained depth, and
    the reference codon's synonymous opportunities. ``min_freq`` is an optional additional
    *codon* frequency floor; position-level allele filtering has already occurred, so it is
    zero by default. Reference stop or ambiguous codons are excluded.
    """
    catalog = Path(catalog)
    counts = read_codon_counts(source, catalog)
    if genome != "all":
        genes = gi.read_gene_info(catalog).lazy().filter(pl.col("genome") == genome).select("gene_id")
        counts = counts.join(genes, on="gene_id", how="semi")
    if min_freq > 0:
        reads = pl.col("count").cast(pl.Float64)
        counts = counts.filter(reads >= reads.sum().over(_KEYS) * min_freq)
    references = pl.scan_parquet(catalog.resolve().with_name("reference_codons.parquet")).with_columns(
        pl.col("reference_codon").replace_strict(_OPPORTUNITIES, default=None, return_dtype=pl.Float64).alias("s_sites")
    ).filter(pl.col("s_sites").is_not_null()).select(*_KEYS, "s_sites")
    packed = (
        counts.with_columns(pl.col("codon").replace_strict(_CODON_IDS, return_dtype=pl.UInt8))
        .group_by(_KEYS)
        .agg(pl.col("codon"), pl.col("count"), pl.col("count").cast(pl.Int128).sum().alias("depth"))
        .filter(pl.col("depth") >= min_cov)
        .with_columns(pl.col("depth").cast(pl.UInt64))
        .join(references, on=_KEYS)
        .sort(_KEYS)
    )
    packed.sink_parquet(destination, engine="streaming", compression="zstd", metadata=metadata or {})


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


def _sum_selection_batches(aligned: pl.LazyFrame, integration: str, batch_size: int) -> pl.LazyFrame:
    """Reduce aligned codons, in order, to one row of selection metrics per gene.

    Codons are added to their gene's running total one at a time, in row order, so the
    floating-point sums do not depend on ``batch_size``.
    """
    totals = {}
    for batch in aligned.collect_batches(chunk_size=batch_size, engine="streaming"):
        if not batch.height:
            continue
        values = _batch_totals(batch, integration)
        gene_ids, index = np.unique(batch["gene_id"].to_numpy(), return_inverse=True)
        running = np.array([totals.get(int(gene), _NO_TOTALS) for gene in gene_ids], dtype=np.float64)
        np.add.at(running, index, values)
        totals.update(zip(map(int, gene_ids), running, strict=True))
    return selection_frame(totals)


def vectorized_pair_selection(left, right, *, allele_integration: str = options.ALLELE_INTEGRATION,
                              batch_size: int = SELECTION_BATCH_ROWS) -> pl.LazyFrame:
    """Compare two packed samples in bounded batches. Results do not depend on ``batch_size``."""
    _selection_arguments(1, 0.0, allele_integration, batch_size)
    a = pl.scan_parquet(left).set_sorted(_KEYS)
    b = pl.scan_parquet(right).drop("s_sites").set_sorted(_KEYS)
    aligned = a.join(b, on=_KEYS, suffix="_2", maintain_order="left")
    return _sum_selection_batches(aligned, allele_integration, batch_size)


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
    reference_file = Path(gene_info).resolve().with_name("reference_codons.parquet")
    with TemporaryDirectory(prefix="zipstrain_reference_selection_") as staging:
        packed = Path(staging)/"sample.parquet"
        _pack_callable_codons(codon_counts, gene_info, packed, min_cov, min_freq=min_freq)
        references = pl.scan_parquet(reference_file).select(*_KEYS, pl.col("reference_codon").replace_strict(
            _CODON_IDS, default=None, return_dtype=pl.UInt8).alias("reference_id"))
        aligned = (
            pl.scan_parquet(packed).set_sorted(_KEYS)
            .rename({"codon": "codon_2", "count": "count_2", "depth": "depth_2"})
            .join(references, on=_KEYS, maintain_order="left")
            .with_columns(pl.concat_list("reference_id").alias("codon"),
                          pl.concat_list(pl.lit(1, dtype=pl.UInt64)).alias("count"),
                          pl.lit(1, dtype=pl.UInt64).alias("depth"))
        )
        return _sum_selection_batches(aligned, allele_integration, batch_size)


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
    with TemporaryDirectory(prefix="zipstrain_pair_selection_") as staging:
        left, right = Path(staging)/"left.parquet", Path(staging)/"right.parquet"
        _pack_callable_codons(codon_counts_1, gene_info, left, min_cov, min_freq=min_freq)
        _pack_callable_codons(codon_counts_2, gene_info, right, min_cov, min_freq=min_freq)
        return vectorized_pair_selection(left, right, allele_integration=allele_integration, batch_size=batch_size)
