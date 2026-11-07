"""Content-addressed staging of pairwise selection for table and matrix comparisons."""

from itertools import combinations
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory

import polars as pl

from zipstrain import codons, gene_info as gi, utils
from zipstrain.resource_limits import map_jobs


_CACHE_VERSION = "matrix-selection-1"
_CACHE_KEY = "zipstrain_matrix_selection_id"
_BYTES_PER_WORST_CASE_ROW = 4096 * 256


def _pack_job(job: tuple) -> None:
    """Pack one sample, atomically, in a worker or in this process."""
    source, catalog, target, min_cov, genome, identity = job
    with TemporaryDirectory(dir=target.parent) as staging:
        temporary = Path(staging) / target.name
        codons._pack_callable_codons(source, catalog, temporary, min_cov, genome, {_CACHE_KEY: identity})
        temporary.replace(target)


def _selection_pair_job(job: tuple) -> None:
    """Calculate and cache one pair's selection result, atomically."""
    left, right, name_a, name_b, allele_integration, batch_size, identity, target = job
    pair = codons.vectorized_pair_selection(
        left, right, allele_integration=allele_integration, batch_size=batch_size
    ).with_columns(sample_1=pl.lit(name_a), sample_2=pl.lit(name_b))
    with TemporaryDirectory(dir=target.parent) as staging:
        temporary = Path(staging) / target.name
        pair.sink_parquet(temporary, engine="streaming", compression="zstd", metadata={_CACHE_KEY: identity})
        temporary.replace(target)


def stage_selection_pairs(profiles, fingerprints, catalog: Path, directory: Path, *,
                          min_cov: int, genome: str, allele_integration: str,
                          memory_limit_gb: float, codon_tables=None, codon_digests=None,
                          pool=None, workers: int = 1) -> list[Path]:
    """Validate samples and cache complete selection results by content.

    Return one result path per sample pair, in sorted-name order. Only samples in pairs
    without a result are packed, and packed samples are removed afterwards, so reruns and
    added samples cost only what is new. ``codon_tables`` maps a sample name to its table when
    it does not sit beside the profile. ``codon_digests`` saves hashing tables the caller
    already hashed. ``pool`` runs the work in ``workers`` processes, which share
    ``memory_limit_gb``. Results do not depend on any of these.
    """
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
        source = (codon_tables or {}).get(name)
        if source is None and not path.name.endswith(suffix):
            raise ValueError(f"Cannot find sibling codon evidence for {path}. Expected a filename ending in {suffix}.")
        source = Path(source) if source is not None else codons.sibling_codon_table(path)
        if not source.is_file():
            raise ValueError(f"Missing observed codon table for {name}. Profile with --collect-codons: {source}")
        if pl.read_parquet_metadata(source).get("zipstrain_source_profile_sha256") != fingerprints[name]:
            raise ValueError(f"Codon evidence does not match its sample position profile: {name}")
        identity = gi.content_digest([
            _CACHE_VERSION, annotation[gi.ANNOTATION_ID], reference_digest, fingerprints[name],
            (codon_digests or {}).get(name) or utils.sha256_file(source), min_cov, genome,
        ])
        samples[name] = source, directory / f"sample-{identity}.parquet", identity

    results, needed = [], []
    for name_a, name_b in combinations(samples, 2):
        identity = gi.content_digest([_CACHE_VERSION, name_a, name_b, samples[name_a][2], samples[name_b][2], allele_integration])
        target = directory / f"pair-{identity}.parquet"
        results.append(target)
        if utils._read_custom_parquet_metadata(target).get(_CACHE_KEY) != identity:
            needed.append((name_a, name_b, identity, target))
    if not needed:
        return results

    # Budget for the worst case of 4,096 observed codon pairs per matched locus,
    # with room for vector indices and intermediate arrays.
    budget = int(memory_limit_gb * 1024**3) // max(1, workers)
    batch_size = max(1, min(codons.SELECTION_BATCH_ROWS, budget // _BYTES_PER_WORST_CASE_ROW))
    used = sorted({name for name_a, name_b, _, _ in needed for name in (name_a, name_b)})
    try:
        map_jobs(pool, _pack_job, [(samples[name][0], catalog, samples[name][1], min_cov, genome, samples[name][2])
                                   for name in used])
        map_jobs(pool, _selection_pair_job, [
            (samples[name_a][1], samples[name_b][1], name_a, name_b, allele_integration, batch_size, identity, target)
            for name_a, name_b, identity, target in needed
        ])
    finally:
        for name in used:
            samples[name][1].unlink(missing_ok=True)
    return results
