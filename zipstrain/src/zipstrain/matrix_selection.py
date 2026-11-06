"""Vectorized selection statistics for the V2 matrix workflow."""

from functools import lru_cache
from itertools import combinations
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import polars as pl

from zipstrain import codons, gene_info as gi, utils


_CODON_IDS = {codon: index for index, codon in enumerate(codons.GENETIC_CODE)}
_BASES = np.array([["ACGT".index(base) for base in codon] for codon in _CODON_IDS], dtype=np.uint8)
_OPPORTUNITIES = {
    codon: codons.synonymous_sites(codon)
    for codon, amino_acid in codons.GENETIC_CODE.items() if amino_acid != "*"
}
_KEYS = ["gene_id", "codon_index"]
_CACHE_VERSION = "matrix-selection-1"
_CACHE_KEY = "zipstrain_matrix_selection_id"


@lru_cache(maxsize=1)
def _consequence_lookup() -> np.ndarray:
    """Return nonsynonymous and synonymous effects by offset and ordered codon pair."""
    effects = np.zeros((3, 64, 64, 2), dtype=np.float64)
    for source in _OPPORTUNITIES:
        for target in _OPPORTUNITIES:
            for offset in range(3):
                if source[offset] != target[offset]:
                    try:
                        synonymous, nonsynonymous = codons.site_changes(source, target, offset)
                    except ValueError:
                        continue
                    effects[offset, _CODON_IDS[source], _CODON_IDS[target]] = nonsynonymous, synonymous
    return effects


def _pack_callable_codons(source: Path, catalog: Path, destination: Path, min_cov: int,
                          genome: str, metadata: dict[str, str]) -> None:
    """Stage one sample as sorted sparse codon groups after validating its evidence."""
    counts = codons.read_codon_counts(source, catalog)
    if genome != "all":
        genes = gi.read_gene_info(catalog).lazy().filter(pl.col("genome") == genome).select("gene_id")
        counts = counts.join(genes, on="gene_id", how="semi")
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
    packed.sink_parquet(destination, engine="streaming", compression="zstd", metadata=metadata)


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


def vectorized_pair_selection(left: Path, right: Path, *, allele_integration: str,
                              batch_size: int = 128) -> pl.LazyFrame:
    """Compare two packed samples in bounded batches using the table-mode estimator."""
    if allele_integration not in codons.ALLELE_INTEGRATION_MODES:
        raise ValueError("allele_integration must be consensus or weighted")
    if batch_size < 1:
        raise ValueError("Selection batch size must be positive")
    a = pl.scan_parquet(left).set_sorted(_KEYS)
    b = pl.scan_parquet(right).drop("s_sites").set_sorted(_KEYS)
    aligned = a.join(b, on=_KEYS, suffix="_2", maintain_order="left")
    totals = {}
    for batch in aligned.collect_batches(chunk_size=batch_size, engine="streaming"):
        if not batch.height:
            continue
        values = _batch_totals(batch, allele_integration)
        gene_ids, index = np.unique(batch["gene_id"].to_numpy(), return_inverse=True)
        grouped = np.zeros((len(gene_ids), 6), dtype=np.float64)
        np.add.at(grouped, index, values)
        for gene, row in zip(gene_ids, grouped, strict=True):
            totals.setdefault(int(gene), np.zeros(6, dtype=np.float64))[:] += row
    return codons.selection_frame(totals)


def stage_selection_pairs(profiles, fingerprints, catalog: Path, directory: Path, *,
                          min_cov: int, genome: str, allele_integration: str,
                          memory_limit_gb: float) -> list[Path]:
    """Validate and cache sparse samples and complete selection results by content."""
    if allele_integration not in codons.ALLELE_INTEGRATION_MODES:
        raise ValueError("allele_integration must be consensus or weighted")
    if min_cov < 1:
        raise ValueError("min_cov must be positive")
    if not isfinite(memory_limit_gb) or memory_limit_gb <= 0:
        raise ValueError("memory_limit_gb must be positive and finite")
    directory.mkdir(parents=True, exist_ok=True)
    annotation = pl.read_parquet_metadata(catalog)
    reference = catalog.resolve().with_name("reference_codons.parquet")
    reference_digest = utils.sha256_file(reference)
    if reference_digest != annotation.get("zipstrain_reference_codons_sha256"):
        raise ValueError("Reference codon contents do not match the preprofile catalog")
    filters = {pl.read_parquet_metadata(path).get("zipstrain_profile_filters") for _, path in profiles}
    if len(filters) > 1:
        raise ValueError("Profiles use incompatible profiling filters")

    samples = {}
    for name, path in sorted(profiles):
        path = Path(path)
        suffix = "_profile.parquet"
        if not path.name.endswith(suffix):
            raise ValueError(f"Cannot find sibling codon evidence for {path}. Expected a filename ending in {suffix}.")
        source = path.with_name(path.name.removesuffix(suffix) + "_codon_counts.parquet")
        if not source.is_file():
            raise ValueError(f"Missing observed codon table for {name}. Profile with --collect-codons: {source}")
        if pl.read_parquet_metadata(source).get("zipstrain_source_profile_sha256") != fingerprints[name]:
            raise ValueError(f"Codon evidence does not match its sample position profile: {name}")
        identity = gi.content_digest([
            _CACHE_VERSION, annotation[gi.ANNOTATION_ID], reference_digest, fingerprints[name],
            utils.sha256_file(source), min_cov, genome,
        ])
        target = directory / f"sample-{identity}.parquet"
        if utils._read_custom_parquet_metadata(target).get(_CACHE_KEY) != identity:
            with TemporaryDirectory(dir=directory) as staging:
                temporary = Path(staging) / target.name
                _pack_callable_codons(source, catalog, temporary, min_cov, genome, {_CACHE_KEY: identity})
                temporary.replace(target)
        samples[name] = (target, identity)

    # Budget for the worst case of 4,096 observed codon pairs per matched locus,
    # with room for vector indices and intermediate arrays.
    batch_size = max(1, min(128, int(memory_limit_gb * 1024**3) // (4096 * 256)))
    results = []
    for name_a, name_b in combinations(samples, 2):
        left, id_a = samples[name_a]
        right, id_b = samples[name_b]
        identity = gi.content_digest([_CACHE_VERSION, name_a, name_b, id_a, id_b, allele_integration])
        target = directory / f"pair-{identity}.parquet"
        if utils._read_custom_parquet_metadata(target).get(_CACHE_KEY) != identity:
            pair = vectorized_pair_selection(
                left, right, allele_integration=allele_integration, batch_size=batch_size
            ).with_columns(sample_1=pl.lit(name_a), sample_2=pl.lit(name_b))
            with TemporaryDirectory(dir=directory) as staging:
                temporary = Path(staging) / target.name
                pair.sink_parquet(temporary, engine="streaming", compression="zstd", metadata={_CACHE_KEY: identity})
                temporary.replace(target)
        results.append(target)
    return results
