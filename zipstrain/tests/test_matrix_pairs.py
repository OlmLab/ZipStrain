from pathlib import Path
from itertools import combinations

from click.testing import CliRunner
import duckdb
import numpy as np
import polars as pl
import pytest

from zipstrain import cli
from zipstrain import compare as cp
from zipstrain import matrix_pairs as mp


def _write_profiles(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    sample_a = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [0, 2, 5],
            "gene": ["gene1", "gene1", "gene2"],
            "genome": ["genome1", "genome1", "genome2"],
            "A": [6, 0, 0],
            "T": [0, 0, 8],
            "C": [0, 3, 0],
            "G": [0, 2, 0],
        }
    )
    sample_b = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chr2"],
            "pos": [1, 5, 7],
            "gene": ["gene1", "gene2", "gene3"],
            "genome": ["genome1", "genome2", "genome2"],
            "A": [0, 1, 0],
            "T": [9, 1, 0],
            "C": [0, 0, 0],
            "G": [0, 0, 6],
        }
    )
    sample_c = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [0, 1, 5, 7],
            "gene": ["gene1", "gene1", "gene2", "gene3"],
            "genome": ["genome1", "genome1", "genome2", "genome2"],
            "A": [5, 0, 0, 0],
            "T": [0, 7, 4, 0],
            "C": [0, 0, 0, 0],
            "G": [0, 0, 0, 5],
        }
    )
    ignored = pl.DataFrame({"cov": [1]})

    sample_a.write_parquet(profile_dir / "sample_a.parquet")
    sample_b.write_parquet(profile_dir / "sample_b.parquet")
    sample_c.write_parquet(profile_dir / "sample_c.parquet")
    ignored.write_parquet(profile_dir / "ignored.parquet")


def _write_profiles_one_based(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    sample_a = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [1, 3, 6],
            "gene": ["gene1", "gene1", "gene2"],
            "genome": ["genome1", "genome1", "genome2"],
            "A": [6, 0, 0],
            "T": [0, 0, 8],
            "C": [0, 3, 0],
            "G": [0, 2, 0],
        }
    )
    sample_b = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2", "chr2"],
            "pos": [2, 6, 8],
            "gene": ["gene1", "gene2", "gene3"],
            "genome": ["genome1", "genome2", "genome2"],
            "A": [0, 1, 0],
            "T": [9, 1, 0],
            "C": [0, 0, 0],
            "G": [0, 0, 6],
        }
    )
    sample_a.write_parquet(profile_dir / "sample_a.parquet")
    sample_b.write_parquet(profile_dir / "sample_b.parquet")


def _write_many_profiles_same_genome(profile_dir: Path, sample_count: int = 5) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    for idx in range(sample_count):
        frame = pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "pos": [0, 1, 2],
                "gene": ["gene1", "gene1", "gene1"],
                "genome": ["genome1", "genome1", "genome1"],
                "A": [6 if idx % 2 == 0 else 0, 0, 6],
                "T": [0, 6 if idx % 2 else 0, 0],
                "C": [0, 0, 0],
                "G": [0, 0, 0],
            }
        )
        frame.write_parquet(profile_dir / f"sample_{idx}.parquet")


def _write_profiles_multiscaffold_same_genome(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    sample_a = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [0, 1, 0, 1],
            "gene": ["gene1", "gene1", "gene2", "gene2"],
            "genome": ["genome1", "genome1", "genome1", "genome1"],
            "A": [6, 6, 6, 6],
            "T": [0, 0, 0, 0],
            "C": [0, 0, 0, 0],
            "G": [0, 0, 0, 0],
        }
    )
    sample_b = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [0, 1, 0, 1],
            "gene": ["gene1", "gene1", "gene2", "gene2"],
            "genome": ["genome1", "genome1", "genome1", "genome1"],
            "A": [7, 7, 7, 7],
            "T": [0, 0, 0, 0],
            "C": [0, 0, 0, 0],
            "G": [0, 0, 0, 0],
        }
    )
    sample_a.write_parquet(profile_dir / "sample_a.parquet")
    sample_b.write_parquet(profile_dir / "sample_b.parquet")


def _write_invalid_append_profile(profile_dir: Path) -> None:
    profile_dir.mkdir(parents=True, exist_ok=True)
    invalid = pl.DataFrame(
        {
            "chrom": ["chr9"],
            "pos": [0],
            "gene": ["gene_bad"],
            "genome": ["genome9"],
            "A": [6],
            "T": [0],
            "C": [0],
            "G": [0],
        }
    )
    invalid.write_parquet(profile_dir / "sample_bad.parquet")


def _write_gene_range_table(gene_range_table: Path) -> None:
    gene_range_table.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(
        {
            "gene": ["gene1", "gene2", "gene3"],
            "scaffold": ["chr1", "chr2", "chr2"],
            "start": [0, 5, 7],
            "end": [2, 5, 7],
        }
    ).write_csv(gene_range_table, separator="\t", include_header=False)


def test_parse_matrix_calculations_supports_all_and_gene_aliases():
    assert mp.parse_matrix_calculations("all") == ("ani", "ibs")
    assert mp.parse_matrix_calculations("all", include_gene_from_all=True) == ("ani", "ibs", "gene")
    assert mp.parse_matrix_calculations("+gene") == ("ani", "gene")
    assert mp.parse_matrix_calculations("popani,max_block") == ("ani", "ibs")


def _load_matrix_compare_db(compare_db: Path):
    conn = duckdb.connect(str(compare_db), read_only=True)
    try:
        metadata = dict(conn.execute("SELECT key, value FROM matrix_compare_metadata ORDER BY key").fetchall())
        completed_pairs = conn.execute(
            """
            SELECT DISTINCT sample_idx_1, sample_idx_2
            FROM matrix_compare_completed_pair_genomes
            ORDER BY sample_idx_1, sample_idx_2
            """
        ).fetchall()
        results = (
            pl.from_arrow(
                conn.execute(
                    """
                    SELECT sample_idx_1, sample_idx_2, sample_1, sample_2, genome_idx, genome,
                           total_positions, share_allele_pos, genome_pop_ani, max_consecutive_length
                    FROM matrix_compare_results
                    ORDER BY sample_idx_1, sample_idx_2, genome_idx
                    """
                ).fetch_arrow_table()
            )
        )
    finally:
        conn.close()
    return metadata, completed_pairs, results


def _load_matrix_compare_gene_results(compare_db: Path) -> pl.DataFrame:
    conn = duckdb.connect(str(compare_db), read_only=True)
    try:
        results = pl.from_arrow(
            conn.execute(
                """
                SELECT sample_idx_1, sample_idx_2, sample_1, sample_2, genome_idx, genome, gene, gene_pop_ani
                FROM matrix_compare_gene_results
                ORDER BY sample_idx_1, sample_idx_2, genome_idx, gene
                """
            ).fetch_arrow_table()
        )
    finally:
        conn.close()
    return results


def _load_matrix_hdf5_store(matrix_hdf5: Path):
    h5py = pytest.importorskip("h5py")
    metadata = mp._load_matrix_hdf5_metadata(matrix_hdf5)
    samples = mp._load_matrix_hdf5_samples(matrix_hdf5)
    genomes = mp._load_matrix_hdf5_genomes(matrix_hdf5)
    scaffolds = mp._load_matrix_hdf5_genome_scaffolds(matrix_hdf5)
    sample_rows = np.arange(len(samples), dtype=np.int64)
    matrix_lengths = {str(spec.genome_idx): spec.matrix_length for spec in genomes}
    numpy_dtype = mp.COUNT_DTYPES[metadata["count_dtype"]]
    with h5py.File(str(matrix_hdf5), "r") as h5_file:
        matrices = {}
        for genome_key, node in h5_file["matrices"].items():
            if hasattr(node, "shape"):
                matrices[genome_key] = np.asarray(node[...])
            else:
                matrices[genome_key] = mp._load_dense_rows_from_sparse_hdf5(
                    indptr_dataset=node["indptr"],
                    indices_dataset=node["indices"],
                    sample_rows=sample_rows,
                    matrix_length=matrix_lengths[genome_key],
                    numpy_dtype=numpy_dtype,
                )
    genes = mp._load_matrix_hdf5_gene_ranges(matrix_hdf5)
    return metadata, samples, genomes, scaffolds, matrices, genes


def _expected_classic_pairwise_gene_results(
    profile_dir: Path,
    *,
    sample_names: tuple[str, ...] = ("sample_a", "sample_b", "sample_c"),
) -> pl.DataFrame:
    expected_frames = []
    for sample_1, sample_2 in combinations(sample_names, 2):
        pair = (
            cp.compare_genes_polars(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                min_gene_compare_len=1,
                genome_scope="all",
                gene_scope="all",
                ani_method="popani",
            )
            .collect(engine="streaming")
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select(["sample_1", "sample_2", "genome", "gene", pl.col("ani").alias("gene_pop_ani")])
        )
        expected_frames.append(pair)
    return pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome", "gene"])


def _expected_classic_pairwise_results(
    profile_dir: Path,
    *,
    calculate: str,
    sample_names: tuple[str, ...] = ("sample_a", "sample_b", "sample_c"),
) -> pl.DataFrame:
    selected_columns = [
        "sample_1",
        "sample_2",
        "genome",
        "total_positions",
        "share_allele_pos",
        "genome_pop_ani",
    ]
    if "ibs" in calculate:
        selected_columns.append("max_consecutive_length")

    expected_frames = []
    for sample_1, sample_2 in combinations(sample_names, 2):
        pair = (
            cp.compare_genomes(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate=calculate,
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select(selected_columns)
        )
        expected_frames.append(pair)
    return pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])


def test_is_hdf5_matrix_store_detects_signature_and_suffix(tmp_path):
    hdf5_path = tmp_path / "matrix.h5"
    hdf5_path.write_bytes(mp.HDF5_FILE_SIGNATURE + b"rest")
    legacy_path = tmp_path / "matrix.duckdb"
    legacy_path.write_bytes(b"DUCK")
    suffix_only = tmp_path / "missing.hdf5"

    assert mp._is_hdf5_matrix_store(hdf5_path) is True
    assert mp._is_hdf5_matrix_store(legacy_path) is False
    assert mp._is_hdf5_matrix_store(suffix_only) is True


def _write_legacy_matrix_db(matrix_db: Path) -> None:
    conn = duckdb.connect(str(matrix_db))
    try:
        mp._init_matrix_db_schema(conn)
        conn.executemany(
            "INSERT INTO matrix_db_metadata VALUES (?, ?)",
            [
                ("profiles_dir", str(matrix_db.parent)),
                ("profile_format", "classic_zipstrain_profile_parquet"),
                ("genome_scope", "all"),
                ("count_dtype", "uint16"),
                ("layout", mp.CURRENT_MATRIX_DB_LAYOUT),
                ("matrix_value_semantics", mp.FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS),
                ("coverage_filter_min_cov", str(mp.MATRIX_BUILD_MIN_COV)),
                ("separator_rows_between_scaffolds", "1"),
            ],
        )
        conn.executemany(
            "INSERT INTO matrix_db_samples VALUES (?, ?, ?)",
            [
                (0, "sample_a", "/tmp/sample_a.parquet"),
                (1, "sample_b", "/tmp/sample_b.parquet"),
                (2, "sample_c", "/tmp/sample_c.parquet"),
            ],
        )
        conn.executemany(
            "INSERT INTO matrix_db_genomes VALUES (?, ?, ?, ?, ?)",
            [
                (0, "genome1", 3, 3, 1),
                (1, "genome2", 3, 3, 1),
            ],
        )
        conn.executemany(
            "INSERT INTO matrix_db_genome_scaffolds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                (0, 0, "genome1", "chr1", 0, 2, 0, 3, 0, 2),
                (1, 0, "genome2", "chr2", 0, 2, 0, 3, 0, 2),
            ],
        )
        matrices = {
            (0, 0): np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1]], dtype=np.uint16),
            (1, 0): np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]], dtype=np.uint16),
            (2, 0): np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]], dtype=np.uint16),
            (0, 1): np.array([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.uint16),
            (1, 1): np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=np.uint16),
            (2, 1): np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=np.uint16),
        }
        for (sample_idx, genome_idx), matrix in matrices.items():
            conn.execute(
                """
                INSERT INTO matrix_db_sample_genome_matrices
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                [
                    sample_idx,
                    genome_idx,
                    "uint16",
                    int(matrix.shape[0]),
                    int(matrix.shape[1]),
                    mp._pack_matrix(matrix),
                ],
            )
    finally:
        conn.close()


def test_export_matrix_db_hdf5_roundtrip(tmp_path):
    h5py = pytest.importorskip("h5py")
    matrix_db = tmp_path / "matrix.duckdb"
    matrix_hdf5 = tmp_path / "matrix.h5"
    _write_legacy_matrix_db(matrix_db)

    exported = mp.export_matrix_db_hdf5(
        matrix_db_file=matrix_db,
        output_file=matrix_hdf5,
    )

    assert exported == matrix_hdf5
    metadata = mp._load_matrix_hdf5_metadata(matrix_hdf5)
    samples = mp._load_matrix_hdf5_samples(matrix_hdf5)
    genomes = mp._load_matrix_hdf5_genomes(matrix_hdf5)
    scaffolds = mp._load_matrix_hdf5_genome_scaffolds(matrix_hdf5)

    assert metadata["input_format"] == "hdf5"
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_LAYOUT
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert [spec.genome for spec in genomes] == ["genome1", "genome2"]
    assert [spec.chrom for spec in scaffolds] == ["chr1", "chr2"]
    with h5py.File(str(matrix_hdf5), "r") as h5_file:
        assert h5_file["samples"]["sample_idx"].maxshape == (None,)
        assert h5_file["samples"]["sample_name"].maxshape == (None,)
        assert h5_file["matrices"]["0"].maxshape == (None, 3, 4)


def test_export_matrix_db_hdf5_sparse_roundtrip(tmp_path):
    h5py = pytest.importorskip("h5py")
    matrix_db = tmp_path / "matrix.duckdb"
    matrix_hdf5 = tmp_path / "matrix_sparse.h5"
    _write_legacy_matrix_db(matrix_db)

    exported = mp.export_matrix_db_hdf5(
        matrix_db_file=matrix_db,
        output_file=matrix_hdf5,
        sparse=True,
    )

    assert exported == matrix_hdf5
    metadata, samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_SPARSE_LAYOUT
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert [spec.genome for spec in genomes] == ["genome1", "genome2"]
    assert [spec.chrom for spec in scaffolds] == ["chr1", "chr2"]
    assert matrices["0"][0].tolist() == [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1]]

    with h5py.File(str(matrix_hdf5), "r") as h5_file:
        assert isinstance(h5_file["matrices"]["0"], h5py.Group)


def test_build_matrix_hdf5(tmp_path):
    h5py = pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    _write_profiles(profile_dir)
    progress_events: list[dict[str, object]] = []

    summary = mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        progress_callback=progress_events.append,
    )

    assert summary.profile_files == 3
    assert summary.sample_count == 3
    assert summary.scaffold_count == 2
    assert summary.count_dtype == "uint16"
    assert summary.stored_rows == 6
    assert progress_events[0]["phase"] == "start"
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 6
    assert progress_events[-1]["total"] == 6
    assert any(event["phase"] == "processing" for event in progress_events)

    metadata = mp._load_matrix_hdf5_metadata(matrix_hdf5)
    samples = mp._load_matrix_hdf5_samples(matrix_hdf5)
    genomes = mp._load_matrix_hdf5_genomes(matrix_hdf5)
    scaffolds = mp._load_matrix_hdf5_genome_scaffolds(matrix_hdf5)

    assert metadata["input_format"] == "hdf5"
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_LAYOUT
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert [spec.genome for spec in genomes] == ["genome1", "genome2"]
    assert [spec.chrom for spec in scaffolds] == ["chr1", "chr2"]

    with h5py.File(str(matrix_hdf5), "r") as h5_file:
        assert h5_file["samples"]["sample_idx"].maxshape == (None,)
        assert h5_file["samples"]["sample_name"].maxshape == (None,)
        assert h5_file["matrices"]["0"].maxshape == (None, 3, 4)
        assert h5_file["matrices"]["1"].maxshape == (None, 3, 4)
        genome1 = np.asarray(h5_file["matrices"]["0"][...])
        genome2 = np.asarray(h5_file["matrices"]["1"][...])

    assert genome1[0].tolist() == [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1]]
    assert genome1[1].tolist() == [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
    assert genome1[2].tolist() == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
    assert genome2[0].tolist() == [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    assert genome2[1].tolist() == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]
    assert genome2[2].tolist() == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]


def test_build_matrix_hdf5_sparse(tmp_path):
    h5py = pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix_sparse.h5"
    _write_profiles(profile_dir)

    summary = mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        sparse=True,
    )

    assert summary.sample_count == 3
    metadata, _samples, _genomes, _scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_SPARSE_LAYOUT
    assert matrices["0"][0].tolist() == [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1]]
    assert matrices["1"][2].tolist() == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]

    with h5py.File(str(matrix_hdf5), "r") as h5_file:
        matrix_node = h5_file["matrices"]["0"]
        assert isinstance(matrix_node, h5py.Group)
        assert matrix_node["indptr"].maxshape == (None,)
        assert matrix_node["indices"].maxshape == (None,)


def test_build_matrix_hdf5_with_gene_ranges(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    _write_profiles(profile_dir)
    _write_gene_range_table(gene_range_table)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )

    metadata, _samples, _genomes, _scaffolds, _matrices, genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["has_gene_ranges"] == "1"
    assert [gene.gene for gene in genes] == ["gene1", "gene2", "gene3"]
    assert [(gene.genome, gene.axis_start, gene.axis_end) for gene in genes] == [
        ("genome1", 0, 2),
        ("genome2", 0, 0),
        ("genome2", 2, 2),
    ]


def test_build_matrix_hdf5_expands_scaffold_span_for_scaffold_relative_gene_ranges(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir(parents=True, exist_ok=True)
    matrix_hdf5 = tmp_path / "matrix.h5"
    gene_range_table = tmp_path / "gene_ranges.tsv"

    pl.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos": [100, 101],
            "gene": ["NA", "NA"],
            "genome": ["genome1", "genome1"],
            "A": [6, 0],
            "T": [0, 6],
            "C": [0, 0],
            "G": [0, 0],
        }
    ).write_parquet(profile_dir / "sample_a.parquet")
    pl.DataFrame(
        {
            "gene": ["gene_early", "gene_observed"],
            "scaffold": ["chr1", "chr1"],
            "start": [10, 100],
            "end": [20, 101],
        }
    ).write_csv(gene_range_table, separator="\t", include_header=False)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )

    metadata, _samples, genomes, scaffolds, matrices, genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["has_gene_ranges"] == "1"
    assert genomes[0].matrix_length == 92
    assert scaffolds[0].min_pos == 10
    assert scaffolds[0].max_pos == 101
    assert [(gene.gene, gene.axis_start, gene.axis_end) for gene in genes] == [
        ("gene_early", 0, 10),
        ("gene_observed", 90, 91),
    ]
    genome_matrix = matrices["0"]
    assert genome_matrix.shape == (1, 92, 4)
    assert genome_matrix[0, 90].tolist() == [1, 0, 0, 0]
    assert genome_matrix[0, 91].tolist() == [0, 1, 0, 0]


def test_build_matrix_hdf5_gene_ranges_follow_multiscaffold_axis_offsets(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    _write_profiles_multiscaffold_same_genome(profile_dir)
    pl.DataFrame(
        {
            "gene": ["gene_chr1", "gene_chr2"],
            "scaffold": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [1, 1],
        }
    ).write_csv(gene_range_table, separator="\t", include_header=False)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )

    metadata, _samples, _genomes, scaffolds, _matrices, genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["separator_rows_between_scaffolds"] == "1"
    scaffold_by_chrom = {spec.chrom: spec for spec in scaffolds}
    gene_by_name = {spec.gene: spec for spec in genes}

    assert gene_by_name["gene_chr1"].axis_start == scaffold_by_chrom["chr1"].axis_start
    assert gene_by_name["gene_chr1"].axis_end == scaffold_by_chrom["chr1"].axis_end
    assert gene_by_name["gene_chr2"].axis_start == scaffold_by_chrom["chr2"].axis_start
    assert gene_by_name["gene_chr2"].axis_end == scaffold_by_chrom["chr2"].axis_end


def test_export_matrix_db_hdf5_reports_progress(tmp_path):
    pytest.importorskip("h5py")
    matrix_db = tmp_path / "matrix.duckdb"
    matrix_hdf5 = tmp_path / "matrix.h5"
    _write_legacy_matrix_db(matrix_db)

    events: list[dict[str, object]] = []
    mp.export_matrix_db_hdf5(
        matrix_db_file=matrix_db,
        output_file=matrix_hdf5,
        progress_callback=events.append,
    )

    assert events[0]["phase"] == "start"
    assert events[-1]["phase"] == "done"
    assert any(event["phase"] == "processing" for event in events)
    advance_events = [event for event in events if event["phase"] == "advance"]
    assert advance_events
    assert advance_events[-1]["stored_rows"] == 6


def test_matrix_compare_rejects_legacy_duckdb_matrix_input(tmp_path):
    profile_dir = tmp_path / "profiles"
    legacy_matrix_db = tmp_path / "matrix.duckdb"
    compare_db = tmp_path / "compare.duckdb"
    _write_profiles(profile_dir)
    _write_legacy_matrix_db(legacy_matrix_db)

    with pytest.raises(ValueError, match="convert it first with 'zipstrain utilities matrix-db-to-hdf5'"):
        mp.matrix_compare(
            matrix_db_file=legacy_matrix_db,
            output_file=compare_db,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            backend="numpy",
            calculate="ani",
        )


def test_matrix_compare_hdf5_gene_results_match_classic_compare(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    compare_db = tmp_path / "compare_hdf5.duckdb"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    _write_profiles(profile_dir)
    _write_gene_range_table(gene_range_table)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )
    summary = mp.matrix_compare(
        matrix_db_file=matrix_hdf5,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="all",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    assert summary.requested_pairs == 3
    metadata, _completed_pairs, _results = _load_matrix_compare_db(compare_db)
    gene_results = _load_matrix_compare_gene_results(compare_db).select(
        ["sample_1", "sample_2", "genome", "gene", "gene_pop_ani"]
    )
    expected = _expected_classic_pairwise_gene_results(profile_dir)

    assert metadata["calculate"] == "ani+ibs+gene"
    assert gene_results.sort(["sample_1", "sample_2", "genome", "gene"]).equals(expected)


def test_matrix_compare_hdf5_gene_results_match_classic_compare_multiscaffold(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    compare_db = tmp_path / "compare_hdf5.duckdb"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    _write_profiles_multiscaffold_same_genome(profile_dir)
    pl.DataFrame(
        {
            "gene": ["gene1", "gene2"],
            "scaffold": ["chr1", "chr2"],
            "start": [0, 0],
            "end": [1, 1],
        }
    ).write_csv(gene_range_table, separator="\t", include_header=False)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )
    mp.matrix_compare(
        matrix_db_file=matrix_hdf5,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="+gene",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    gene_results = _load_matrix_compare_gene_results(compare_db).select(
        ["sample_1", "sample_2", "genome", "gene", "gene_pop_ani"]
    )
    expected = _expected_classic_pairwise_gene_results(
        profile_dir,
        sample_names=("sample_a", "sample_b"),
    )

    assert gene_results.sort(["sample_1", "sample_2", "genome", "gene"]).equals(expected)


def test_accumulate_gene_counts_from_full_torch_masks_respects_gene_boundaries_with_gaps():
    torch = pytest.importorskip("torch")
    gene_ranges = [
        mp.GeneRangeSpec(
            gene_idx=0,
            genome_idx=0,
            genome="g",
            chrom="chr1",
            gene="gene1",
            axis_start=0,
            axis_end=1,
        ),
        mp.GeneRangeSpec(
            gene_idx=1,
            genome_idx=0,
            genome="g",
            chrom="chr2",
            gene="gene2",
            axis_start=3,
            axis_end=4,
        ),
    ]
    total_mask = torch.tensor(
        [
            [1, 1],
            [1, 0],
            [0, 0],
            [1, 1],
            [0, 1],
        ],
        dtype=torch.bool,
    )
    shared_mask = torch.tensor(
        [
            [1, 0],
            [0, 0],
            [0, 0],
            [1, 1],
            [0, 1],
        ],
        dtype=torch.bool,
    )

    gene_total, gene_shared = mp._accumulate_gene_counts_from_full_torch_masks(
        torch_module=torch,
        total_mask=total_mask,
        shared_mask=shared_mask,
        gene_ranges=gene_ranges,
    )

    np.testing.assert_array_equal(
        gene_total,
        np.array(
            [
                [2, 1],
                [1, 2],
            ],
            dtype=np.int64,
        ),
    )
    np.testing.assert_array_equal(
        gene_shared,
        np.array(
            [
                [1, 0],
                [1, 2],
            ],
            dtype=np.int64,
        ),
    )


def test_max_ibs_from_shared_mask_numpy_matches_numpy_reference():
    shared_mask = np.array(
        [
            [True, False, True],
            [True, True, True],
            [False, True, False],
            [True, True, True],
            [True, False, True],
        ],
        dtype=bool,
    )
    _current_runs, expected = mp._update_ibs_numpy(
        shared_mask=shared_mask,
        current_runs=np.zeros(shared_mask.shape[1], dtype=np.int64),
        max_runs=np.zeros(shared_mask.shape[1], dtype=np.int64),
    )

    observed = mp._max_ibs_from_shared_mask_numpy(shared_mask)

    np.testing.assert_array_equal(observed, expected)


def test_matrix_compare_hdf5_explicit_gene_results_match_classic_compare(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    compare_db = tmp_path / "compare_hdf5_gene.duckdb"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    _write_profiles(profile_dir)
    _write_gene_range_table(gene_range_table)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )
    summary = mp.matrix_compare(
        matrix_db_file=matrix_hdf5,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="+gene",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    assert summary.requested_pairs == 3
    metadata, _completed_pairs, results = _load_matrix_compare_db(compare_db)
    gene_results = _load_matrix_compare_gene_results(compare_db).select(
        ["sample_1", "sample_2", "genome", "gene", "gene_pop_ani"]
    )
    expected = _expected_classic_pairwise_gene_results(profile_dir)

    assert metadata["calculate"] == "ani+gene"
    assert results.select(["sample_1", "sample_2", "genome"]).height == 3
    assert gene_results.sort(["sample_1", "sample_2", "genome", "gene"]).equals(expected)


def test_matrix_compare_gene_requires_gene_annotations(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    compare_db = tmp_path / "compare_hdf5_gene.duckdb"
    _write_profiles(profile_dir)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )

    with pytest.raises(ValueError, match="does not contain gene annotations"):
        mp.matrix_compare(
            matrix_db_file=matrix_hdf5,
            output_file=compare_db,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            backend="torch-cpu",
            calculate="+gene",
            loader_executor_kind="thread",
            writer_executor_kind="thread",
        )


def test_matrix_compare_torch_cpu_matches_classic_compare(tmp_path):
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare_torch.duckdb"
    _write_profiles(profile_dir)

    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="ani",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    assert summary.requested_pairs == 3
    _metadata, _completed_pairs, results = _load_matrix_compare_db(output_file)
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])
    expected = _expected_classic_pairwise_results(profile_dir, calculate="ani")

    assert actual.equals(expected)


def test_matrix_compare_torch_cpu_with_ibs_matches_classic_compare(tmp_path):
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare_torch_ibs.duckdb"
    _write_profiles(profile_dir)

    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="ani+ibs",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    assert summary.requested_pairs == 3
    _metadata, _completed_pairs, results = _load_matrix_compare_db(output_file)
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani", "max_consecutive_length"]
    ).sort(["sample_1", "sample_2", "genome"])
    expected = _expected_classic_pairwise_results(profile_dir, calculate="ani+ibs")

    assert actual.equals(expected)


def test_matrix_compare_direct_hdf5_torch_cpu_matches_classic_compare(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare_hdf5_classic.duckdb"
    _write_profiles(profile_dir)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_hdf5,
        count_dtype="uint16",
        memory_limit_gb=1.0,
    )
    summary = mp.matrix_compare(
        matrix_db_file=matrix_hdf5,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="ani+ibs",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    assert summary.requested_pairs == 3
    _metadata, _completed_pairs, results = _load_matrix_compare_db(output_file)
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani", "max_consecutive_length"]
    ).sort(["sample_1", "sample_2", "genome"])
    expected = _expected_classic_pairwise_results(profile_dir, calculate="ani+ibs")

    assert actual.equals(expected)


def test_build_matrix_db(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    _write_profiles(profile_dir)
    progress_events: list[dict[str, object]] = []

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        progress_callback=progress_events.append,
    )

    assert summary.profile_files == 3
    assert summary.sample_count == 3
    assert summary.scaffold_count == 2
    assert summary.count_dtype == "uint16"
    assert summary.stored_rows == 6
    assert progress_events[0]["phase"] == "start"
    assert any(event["phase"] == "start" for event in progress_events)
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 6
    assert progress_events[-1]["total"] == 6
    assert any(event["phase"] == "processing" for event in progress_events)
    assert any(event["sample_name"] == "sample_a" for event in progress_events)
    assert all(event["scaffold"] == "" for event in progress_events if "scaffold" in event)

    metadata, samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert [(spec.genome_idx, spec.genome, spec.matrix_length, spec.true_length, spec.scaffold_count) for spec in genomes] == [
        (0, "genome1", 3, 3, 1),
        (1, "genome2", 3, 3, 1),
    ]
    assert [(spec.genome_idx, spec.scaffold_ordinal, spec.genome, spec.chrom, spec.axis_start, spec.axis_end, spec.vector_length) for spec in scaffolds] == [
        (0, 0, "genome1", "chr1", 0, 2, 3),
        (1, 0, "genome2", "chr2", 0, 2, 3),
    ]
    assert metadata["matrix_value_semantics"] == mp.FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS
    assert metadata["coverage_filter_min_cov"] == str(mp.MATRIX_BUILD_MIN_COV)
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_LAYOUT
    assert metadata["separator_rows_between_scaffolds"] == "1"

    assert matrices["0"][0].tolist() == [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1]]
    assert matrices["0"][1].tolist() == [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
    assert matrices["0"][2].tolist() == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
    assert matrices["1"][0].tolist() == [[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    assert matrices["1"][1].tolist() == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]
    assert matrices["1"][2].tolist() == [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]


def test_build_matrix_db_with_small_commit_batches(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_small_batches.h5"
    _write_profiles(profile_dir)

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        count_dtype="uint16",
        memory_limit_gb=1.0,
        commit_batch_gb=1e-9,
    )

    assert summary.stored_rows == 6
    _metadata, samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_db)
    assert len(samples) == 3
    assert len(genomes) == 2
    assert len(scaffolds) == 2
    assert matrices["0"].shape[0] == 3
    assert matrices["1"].shape[0] == 3


def test_build_matrix_db_with_optional_bed_file(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_bed.h5"
    bed_file = tmp_path / "genomes.bed"
    _write_profiles_one_based(profile_dir)
    bed_file.write_text("chr1\t0\t2\nchr1\t2\t5\nchr2\t5\t7\nchr2\t7\t10\n")

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        bed_file=bed_file,
        count_dtype="uint16",
        memory_limit_gb=1.0,
    )

    assert summary.profile_files == 2
    assert summary.scaffold_count == 2

    _metadata, _samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_db)
    assert [(spec.genome_idx, spec.genome, spec.matrix_length, spec.true_length, spec.scaffold_count) for spec in genomes] == [
        (0, "genome1", 5, 5, 1),
        (1, "genome2", 5, 5, 1),
    ]
    assert [(spec.genome_idx, spec.scaffold_ordinal, spec.genome, spec.chrom, spec.axis_start, spec.axis_end, spec.vector_length) for spec in scaffolds] == [
        (0, 0, "genome1", "chr1", 0, 4, 5),
        (1, 0, "genome2", "chr2", 0, 4, 5),
    ]
    assert matrices["0"][0].tolist() == [
        [1, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 1, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
    assert matrices["0"][1].tolist() == [
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]


def test_append_matrix_db_success(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_db = tmp_path / "matrix.h5"

    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    initial_summary = mp.build_matrix_db(
        profile_dir=initial_profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    assert initial_summary.sample_count == 2

    progress_events: list[dict[str, object]] = []
    append_summary = mp.append_matrix_db(
        profile_dir=append_profile_dir,
        matrix_db_file=matrix_db,
        memory_limit_gb=1.0,
        progress_callback=progress_events.append,
    )

    assert append_summary.appended_sample_count == 1
    assert append_summary.total_sample_count == 3
    assert append_summary.stored_rows == 2
    assert progress_events[0]["phase"] == "start"
    assert any(event["phase"] == "start" for event in progress_events)
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 2
    assert progress_events[-1]["total"] == 2

    _metadata, samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert len(genomes) == 2
    assert len(scaffolds) == 2
    assert matrices["0"].shape[0] == 3
    assert matrices["1"].shape[0] == 3


def test_append_matrix_hdf5_success(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_hdf5 = tmp_path / "matrix.h5"

    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    initial_summary = mp.build_matrix_hdf5(
        profile_dir=initial_profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )
    assert initial_summary.sample_count == 2

    progress_events: list[dict[str, object]] = []
    append_summary = mp.append_matrix_hdf5(
        profile_dir=append_profile_dir,
        matrix_hdf5_file=matrix_hdf5,
        memory_limit_gb=1.0,
        progress_callback=progress_events.append,
    )

    assert append_summary.appended_sample_count == 1
    assert append_summary.total_sample_count == 3
    assert append_summary.stored_rows == 2
    assert progress_events[0]["phase"] == "start"
    assert progress_events[-1]["phase"] == "done"
    assert progress_events[-1]["completed"] == 2
    assert progress_events[-1]["total"] == 2

    metadata, samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["input_format"] == "hdf5"
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert len(genomes) == 2
    assert len(scaffolds) == 2
    assert matrices["0"].shape[0] == 3
    assert matrices["1"].shape[0] == 3


def test_append_matrix_hdf5_sparse_success(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_hdf5 = tmp_path / "matrix_sparse.h5"

    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    initial_summary = mp.build_matrix_hdf5(
        profile_dir=initial_profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
        sparse=True,
    )
    assert initial_summary.sample_count == 2

    append_summary = mp.append_matrix_hdf5(
        profile_dir=append_profile_dir,
        matrix_hdf5_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )

    assert append_summary.appended_sample_count == 1
    metadata, samples, _genomes, _scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_SPARSE_LAYOUT
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert matrices["0"].shape[0] == 3
    assert matrices["1"].shape[0] == 3


def test_append_matrix_hdf5_uses_in_place_path_for_resizable_store(tmp_path, monkeypatch):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_hdf5 = tmp_path / "matrix.h5"

    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    mp.build_matrix_hdf5(
        profile_dir=initial_profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )

    def fail_rewrite(**kwargs):
        raise AssertionError("rewrite fallback should not be used for a resizable HDF5 store")

    monkeypatch.setattr(mp, "_append_matrix_hdf5_via_rewrite", fail_rewrite)
    summary = mp.append_matrix_hdf5(
        profile_dir=append_profile_dir,
        matrix_hdf5_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )

    assert summary.appended_sample_count == 1


def test_append_matrix_db_rejects_incompatible_profile_without_mutation(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_bad"
    matrix_db = tmp_path / "matrix.h5"

    _write_profiles(initial_profile_dir)
    (initial_profile_dir / "sample_c.parquet").unlink()
    _write_invalid_append_profile(append_profile_dir)

    mp.build_matrix_db(
        profile_dir=initial_profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    metadata_before, samples_before, genomes_before, scaffolds_before, matrices_before, genes_before = _load_matrix_hdf5_store(matrix_db)

    with pytest.raises(ValueError, match="does not match the existing matrix store contract"):
        mp.append_matrix_db(
            profile_dir=append_profile_dir,
            matrix_db_file=matrix_db,
            memory_limit_gb=1.0,
        )

    metadata_after, samples_after, genomes_after, scaffolds_after, matrices_after, genes_after = _load_matrix_hdf5_store(matrix_db)
    assert metadata_after == metadata_before
    assert samples_after == samples_before
    assert genomes_after == genomes_before
    assert scaffolds_after == scaffolds_before
    assert matrices_after.keys() == matrices_before.keys()
    for genome_idx in matrices_before:
        assert np.array_equal(matrices_after[genome_idx], matrices_before[genome_idx])
    assert genes_after == genes_before


def test_append_matrix_hdf5_rejects_incompatible_profile_without_mutation(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_bad"
    matrix_hdf5 = tmp_path / "matrix.h5"

    _write_profiles(initial_profile_dir)
    (initial_profile_dir / "sample_c.parquet").unlink()
    _write_invalid_append_profile(append_profile_dir)

    mp.build_matrix_hdf5(
        profile_dir=initial_profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )
    metadata_before, samples_before, genomes_before, scaffolds_before, matrices_before, _genes_before = _load_matrix_hdf5_store(matrix_hdf5)

    with pytest.raises(ValueError, match="does not match the existing matrix store contract"):
        mp.append_matrix_hdf5(
            profile_dir=append_profile_dir,
            matrix_hdf5_file=matrix_hdf5,
            memory_limit_gb=1.0,
        )

    metadata_after, samples_after, genomes_after, scaffolds_after, matrices_after, _genes_after = _load_matrix_hdf5_store(matrix_hdf5)
    assert metadata_after == metadata_before
    assert samples_after == samples_before
    assert genomes_after == genomes_before
    assert scaffolds_after == scaffolds_before
    assert matrices_after.keys() == matrices_before.keys()
    for genome_key in matrices_before:
        assert np.array_equal(matrices_after[genome_key], matrices_before[genome_key])


def test_matrix_compare_matches_pairwise_compare(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    progress_events: list[dict[str, object]] = []

    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
        progress_callback=progress_events.append,
    )

    assert summary.requested_pairs == 3
    assert summary.backend == "numpy"
    assert summary.anchor_groups == 2
    assert progress_events
    assert all(event["phase"] == "advance" for event in progress_events)
    assert progress_events[-1]["completed"] == 6
    assert progress_events[-1]["total"] == 6
    assert any(event["anchor_name"] == "sample_a" for event in progress_events)
    assert all(event["scaffold"] == "" for event in progress_events if "scaffold" in event)
    metadata, _completed_pairs, results = _load_matrix_compare_db(output_file)
    assert metadata["calculate"] == "ani"
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select([
                "sample_1",
                "sample_2",
                "genome",
                "total_positions",
                "share_allele_pos",
                "genome_pop_ani",
            ])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])

    assert actual.equals(expected)


def test_matrix_compare_with_ibs_matches_pairwise_compare(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare_ibs.duckdb"
    _write_profiles(profile_dir)

    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani+ibs",
    )

    assert summary.requested_pairs == 3
    metadata, _completed_pairs, results = _load_matrix_compare_db(output_file)
    assert metadata["calculate"] == "ani+ibs"
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani", "max_consecutive_length"]
    ).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani+ibs",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(
                sample_1=pl.lit(sample_1),
                sample_2=pl.lit(sample_2),
            )
            .select([
                "sample_1",
                "sample_2",
                "genome",
                "total_positions",
                "share_allele_pos",
                "genome_pop_ani",
                "max_consecutive_length",
            ])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])

    assert actual.equals(expected)


def test_matrix_compare_loads_targets_in_batches(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    original_load_indices = mp._Hdf5GenomeMatrixNumpyDataset.load_indices
    call_sizes: list[int] = []

    def tracking_load_indices(self, sample_rows):
        call_sizes.append(len(sample_rows))
        return original_load_indices(self, sample_rows)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    monkeypatch.setattr(mp._Hdf5GenomeMatrixNumpyDataset, "load_indices", tracking_load_indices)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)

    mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )

    assert call_sizes
    assert sorted(call_sizes) == [1, 1, 2, 2]


def test_cli_matrix_build_and_compare(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    export_file = tmp_path / "matrix_compare.parquet"
    _write_profiles(profile_dir)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
            "--export-batch-mb",
            "64",
        ],
    )
    assert build_result.exit_code == 0
    assert "scaffolds=2" in build_result.output

    compare_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare",
            "--matrix-db-file",
            str(matrix_db),
            "--output-file",
            str(output_file),
            "--memory-limit-gb",
            "1",
            "--target-queue-size",
            "1",
            "--loader-executor",
            "thread",
            "--writer-executor",
            "thread",
            "--backend",
            "torch-cpu",
            "--calculate",
            "ani",
        ],
    )
    assert compare_result.exit_code == 0
    assert output_file.exists()
    assert "requested_pairs=3" in compare_result.output
    assert "anchor_groups=2" in compare_result.output
    assert "MATRIX-COMPARE PROGRESS" in compare_result.output
    assert "elapsed=" in compare_result.output
    assert "batch_pairs=" in compare_result.output
    assert "batch_rows=" in compare_result.output

    export_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare-export",
            "--matrix-compare-db-file",
            str(output_file),
            "--output-file",
            str(export_file),
        ],
    )
    assert export_result.exit_code == 0
    out = pl.read_parquet(export_file)
    assert "max_consecutive_length" not in out.columns

    missing_gene_export = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare-export",
            "--matrix-compare-db-file",
            str(output_file),
            "--output-file",
            str(tmp_path / "matrix_compare_gene.parquet"),
            "--table",
            "gene",
        ],
    )
    assert missing_gene_export.exit_code != 0
    assert "does not contain any gene comparison rows" in missing_gene_export.output


def test_cli_matrix_build_sparse_and_compare(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_sparse.h5"
    output_file = tmp_path / "matrix_compare_sparse.duckdb"
    _write_profiles(profile_dir)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
            "--sparse",
        ],
    )
    assert build_result.exit_code == 0

    metadata = mp._load_matrix_hdf5_metadata(matrix_db)
    assert metadata["layout"] == mp.CURRENT_MATRIX_HDF5_SPARSE_LAYOUT

    compare_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare",
            "--matrix-db-file",
            str(matrix_db),
            "--output-file",
            str(output_file),
            "--memory-limit-gb",
            "1",
            "--target-queue-size",
            "1",
            "--loader-executor",
            "thread",
            "--writer-executor",
            "thread",
            "--backend",
            "torch-cpu",
            "--calculate",
            "ani",
        ],
    )
    assert compare_result.exit_code == 0
    assert output_file.exists()


def test_cli_matrix_compare_export_gene_table(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    export_file = tmp_path / "matrix_compare_gene.parquet"
    _write_profiles(profile_dir)
    _write_gene_range_table(gene_range_table)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
            "--gene-range-table",
            str(gene_range_table),
        ],
    )
    assert build_result.exit_code == 0

    compare_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare",
            "--matrix-db-file",
            str(matrix_db),
            "--output-file",
            str(output_file),
            "--memory-limit-gb",
            "1",
            "--target-queue-size",
            "1",
            "--loader-executor",
            "thread",
            "--writer-executor",
            "thread",
            "--backend",
            "torch-cpu",
            "--calculate",
            "+gene",
        ],
    )
    assert compare_result.exit_code == 0

    export_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "matrix-compare-export",
            "--matrix-compare-db-file",
            str(output_file),
            "--output-file",
            str(export_file),
            "--table",
            "gene",
        ],
    )
    assert export_result.exit_code == 0
    out = pl.read_parquet(export_file)
    assert out.columns == [
        "sample_1",
        "sample_2",
        "genome",
        "gene",
        "gene_pop_ani",
    ]
    assert "gene1" in out["gene"].to_list()


def test_export_matrix_compare_parquet_supports_genome_and_gene_tables(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    gene_range_table = tmp_path / "gene_ranges.tsv"
    matrix_db = tmp_path / "matrix.h5"
    compare_db = tmp_path / "matrix_compare.duckdb"
    genome_export = tmp_path / "matrix_compare.parquet"
    gene_export = tmp_path / "matrix_compare_gene.parquet"
    _write_profiles(profile_dir)
    _write_gene_range_table(gene_range_table)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
        gene_range_table=gene_range_table,
    )
    mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="all",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    written_genome = mp.export_matrix_compare_parquet(
        matrix_compare_db_file=compare_db,
        output_file=genome_export,
        table="genome",
    )
    written_gene = mp.export_matrix_compare_parquet(
        matrix_compare_db_file=compare_db,
        output_file=gene_export,
        table="gene",
    )

    assert written_genome == genome_export.resolve()
    assert written_gene == gene_export.resolve()
    genome_df = pl.read_parquet(genome_export)
    gene_df = pl.read_parquet(gene_export)
    assert genome_df.columns[:3] == ["sample_1", "sample_2", "genome"]
    assert gene_df.columns == ["sample_1", "sample_2", "genome", "gene", "gene_pop_ani"]


def test_export_matrix_compare_parquet_gene_requires_gene_rows(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    compare_db = tmp_path / "matrix_compare.duckdb"
    gene_export = tmp_path / "matrix_compare_gene.parquet"
    _write_profiles(profile_dir)

    mp.build_matrix_hdf5(
        profile_dir=profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="ani",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )

    with pytest.raises(ValueError, match="does not contain any gene comparison rows"):
        mp.export_matrix_compare_parquet(
            matrix_compare_db_file=compare_db,
            output_file=gene_export,
            table="gene",
        )


def test_cli_legacy_build_matrix_hdf5_command_is_removed(tmp_path):
    pytest.importorskip("h5py")
    profile_dir = tmp_path / "profiles"
    matrix_hdf5 = tmp_path / "matrix.h5"
    _write_profiles(profile_dir)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-hdf5",
            "--profile-dir",
            str(profile_dir),
            "--output-file",
            str(matrix_hdf5),
        ],
    )
    assert build_result.exit_code != 0
    assert "No such command 'build-matrix-hdf5'" in build_result.output


def test_cli_append_matrix_db(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_db = tmp_path / "matrix.h5"
    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(initial_profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
            "--export-batch-mb",
            "64",
        ],
    )
    assert build_result.exit_code == 0

    append_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "append-matrix-db",
            "--profile-dir",
            str(append_profile_dir),
            "--matrix-db-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert append_result.exit_code == 0
    assert "appended_samples=1" in append_result.output
    assert "total_samples=3" in append_result.output

    _metadata, samples, _genomes, _scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b", "sample_c"]
    assert matrices["0"].shape[0] == 3


def test_cli_legacy_append_matrix_hdf5_command_is_removed(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_hdf5 = tmp_path / "matrix.h5"
    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(initial_profile_dir),
            "--output-file",
            str(matrix_hdf5),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert build_result.exit_code == 0

    append_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "append-matrix-hdf5",
            "--profile-dir",
            str(append_profile_dir),
            "--matrix-hdf5-file",
            str(matrix_hdf5),
        ],
    )
    assert append_result.exit_code != 0
    assert "No such command 'append-matrix-hdf5'" in append_result.output


def test_matrix_compare_resumes_after_matrix_append(tmp_path):
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_db = tmp_path / "matrix.duckdb"
    compare_db = tmp_path / "matrix_compare.duckdb"
    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    mp.build_matrix_db(
        profile_dir=initial_profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )
    first_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )
    assert first_summary.requested_pairs == 1

    mp.append_matrix_db(
        profile_dir=append_profile_dir,
        matrix_db_file=matrix_db,
        memory_limit_gb=1.0,
    )
    second_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )
    assert second_summary.requested_pairs == 2

    _metadata, completed_pairs, results = _load_matrix_compare_db(compare_db)
    assert completed_pairs == [(0, 1), (0, 2), (1, 2)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=(append_profile_dir if sample_1 == "sample_c" else initial_profile_dir) / f"{sample_1}.parquet",
                mpile_contig_2=(append_profile_dir if sample_2 == "sample_c" else initial_profile_dir) / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(sample_1=pl.lit(sample_1), sample_2=pl.lit(sample_2))
            .select(["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)


def test_matrix_compare_resumes_after_matrix_hdf5_append(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("torch")
    initial_profile_dir = tmp_path / "profiles_initial"
    append_profile_dir = tmp_path / "profiles_append"
    matrix_hdf5 = tmp_path / "matrix.h5"
    compare_db = tmp_path / "matrix_compare.duckdb"
    _write_profiles(initial_profile_dir)
    append_profile_dir.mkdir(parents=True, exist_ok=True)
    pl.read_parquet(initial_profile_dir / "sample_c.parquet").write_parquet(append_profile_dir / "sample_c.parquet")
    (initial_profile_dir / "sample_c.parquet").unlink()

    mp.build_matrix_hdf5(
        profile_dir=initial_profile_dir,
        output_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )
    first_summary = mp.matrix_compare(
        matrix_db_file=matrix_hdf5,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="ani",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )
    assert first_summary.requested_pairs == 1

    mp.append_matrix_hdf5(
        profile_dir=append_profile_dir,
        matrix_hdf5_file=matrix_hdf5,
        memory_limit_gb=1.0,
    )
    second_summary = mp.matrix_compare(
        matrix_db_file=matrix_hdf5,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch-cpu",
        calculate="ani",
        loader_executor_kind="thread",
        writer_executor_kind="thread",
    )
    assert second_summary.requested_pairs == 2

    _metadata, completed_pairs, results = _load_matrix_compare_db(compare_db)
    assert completed_pairs == [(0, 1), (0, 2), (1, 2)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=(append_profile_dir if sample_1 == "sample_c" else initial_profile_dir) / f"{sample_1}.parquet",
                mpile_contig_2=(append_profile_dir if sample_2 == "sample_c" else initial_profile_dir) / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(sample_1=pl.lit(sample_1), sample_2=pl.lit(sample_2))
            .select(["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)


def test_matrix_compare_skips_work_when_everything_is_already_done(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    compare_db = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    first_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )
    second_summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=compare_db,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )

    assert first_summary.requested_pairs == 3
    assert second_summary.requested_pairs == 0
    _metadata, completed_pairs, results = _load_matrix_compare_db(compare_db)
    assert completed_pairs == [(0, 1), (0, 2), (1, 2)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"]
    ).sort(["sample_1", "sample_2", "genome"])

    expected_frames = []
    for sample_1, sample_2 in [("sample_a", "sample_b"), ("sample_a", "sample_c"), ("sample_b", "sample_c")]:
        pair = (
            cp.compare_genomes(
                mpile_contig_1=profile_dir / f"{sample_1}.parquet",
                mpile_contig_2=profile_dir / f"{sample_2}.parquet",
                min_cov=mp.MATRIX_BUILD_MIN_COV,
                genome_scope="all",
                ani_method="popani",
                engine="polars",
                calculate="ani",
                stb_file=None,
            )
            .collect(engine="streaming")
            .with_columns(sample_1=pl.lit(sample_1), sample_2=pl.lit(sample_2))
            .select(["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani"])
        )
        expected_frames.append(pair)
    expected = pl.concat(expected_frames).sort(["sample_1", "sample_2", "genome"])
    assert actual.equals(expected)


def test_cli_append_matrix_db_failure_keeps_existing_db(tmp_path):
    pytest.importorskip("h5py")
    initial_profile_dir = tmp_path / "profiles_initial"
    bad_profile_dir = tmp_path / "profiles_bad"
    matrix_db = tmp_path / "matrix.h5"
    _write_profiles(initial_profile_dir)
    (initial_profile_dir / "sample_c.parquet").unlink()
    _write_invalid_append_profile(bad_profile_dir)

    runner = CliRunner()
    build_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "build-matrix-db",
            "--profile-dir",
            str(initial_profile_dir),
            "--output-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
            "--export-batch-mb",
            "64",
        ],
    )
    assert build_result.exit_code == 0
    metadata_before, samples_before, genomes_before, scaffolds_before, matrices_before, genes_before = _load_matrix_hdf5_store(matrix_db)

    append_result = runner.invoke(
        cli.cli,
        [
            "utilities",
            "append-matrix-db",
            "--profile-dir",
            str(bad_profile_dir),
            "--matrix-db-file",
            str(matrix_db),
            "--memory-limit-gb",
            "1",
        ],
    )
    assert append_result.exit_code != 0
    assert append_result.exception is not None
    assert "does not match the existing matrix store contract" in str(append_result.exception)

    metadata_after, samples_after, genomes_after, scaffolds_after, matrices_after, genes_after = _load_matrix_hdf5_store(matrix_db)
    assert metadata_after == metadata_before
    assert samples_after == samples_before
    assert genomes_after == genomes_before
    assert scaffolds_after == scaffolds_before
    assert genes_after == genes_before
    assert matrices_after.keys() == matrices_before.keys()
    for genome_idx in matrices_before:
        assert np.array_equal(matrices_after[genome_idx], matrices_before[genome_idx])


def test_matrix_compare_torch_backend_requires_torch(tmp_path, monkeypatch):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    def _raise_import_error(_name):
        raise ImportError("torch missing")

    monkeypatch.setattr(mp.importlib, "import_module", _raise_import_error)

    with pytest.raises(mp.TorchBackendMissingError, match="Torch backend requested"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            backend="torch",
            memory_limit_gb=1.0,
            calculate="ani",
        )


def test_matrix_compare_ignores_requested_min_cov(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=1,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani",
    )

    assert summary.requested_pairs == 3
    _metadata, _completed_pairs, _results = _load_matrix_compare_db(output_file)




def test_matrix_compare_rejects_unknown_io_executor_kind(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.duckdb"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_profiles(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    with pytest.raises(ValueError, match="loader_executor_kind must be one of"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            loader_executor_kind="bad",
            backend="numpy",
            calculate="ani",
        )

    with pytest.raises(ValueError, match="writer_executor_kind must be one of"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            writer_executor_kind="bad",
            backend="numpy",
            calculate="ani",
        )

    with pytest.raises(ValueError, match="result_transfer_batch_size must be >= 1"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            result_transfer_batch_size=0,
            backend="numpy",
            calculate="ani",
        )


def test_matrix_compare_torch_reuses_target_chunks_across_anchors(tmp_path, monkeypatch):
    torch_module = pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = torch_module

    original_target_load = mp._load_target_queue_block_for_hdf5_torch
    block_sizes: list[int] = []

    def tracking_target_load(*args, **kwargs):
        block_rows = kwargs.get("block_rows")
        if block_rows is None:
            block_rows = args[2]
        block_sizes.append(len(block_rows))
        return original_target_load(*args, **kwargs)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    def identity_prepare(compute_backend, matrix, matrix_value_semantics):
        return matrix

    def zero_compare(
        compute_backend,
        anchor_torch,
        target_torch,
        vector_length: int,
        matrix_value_semantics: str,
        need_ibs: bool = False,
        gene_ranges=None,
    ) -> tuple[np.ndarray, np.ndarray, None, None, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
            None,
            None,
        )

    def download_stub(compute_backend, totals_tensors, shared_tensors):
        max_len = max((len(tensor) for tensor in totals_tensors), default=0)
        out = np.zeros((len(totals_tensors), 2, max_len), dtype=np.int64)
        for idx, (totals_tensor, shared_tensor) in enumerate(zip(totals_tensors, shared_tensors)):
            out[idx, 0, : len(totals_tensor)] = totals_tensor
            out[idx, 1, : len(shared_tensor)] = shared_tensor
        return out

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "_load_target_queue_block_for_hdf5_torch", tracking_target_load)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)
    monkeypatch.setattr(mp, "_prepare_torch_matrix", identity_prepare)
    monkeypatch.setattr(mp, "_compare_anchor_against_target_chunk_torch_device", zero_compare)
    monkeypatch.setattr(mp, "_download_torch_result_tensor_batch", download_stub)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch",
        calculate="ani",
    )

    assert summary.target_chunks == 3
    assert sorted(block_sizes) == [1, 2, 2]


def test_matrix_compare_torch_anchor_queue_batches_host_loads(tmp_path, monkeypatch):
    torch_module = pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = torch_module

    original_anchor_load = mp._load_anchor_queue_batch_for_hdf5_torch
    batch_sizes: list[int] = []

    def tracking_anchor_load(*args, **kwargs):
        batch_rows = kwargs.get("batch_rows")
        if batch_rows is None:
            batch_rows = args[2]
        batch_sizes.append(len(batch_rows))
        return original_anchor_load(*args, **kwargs)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    def identity_prepare(compute_backend, matrix, matrix_value_semantics):
        return matrix

    def zero_compare(
        compute_backend,
        anchor_torch,
        target_torch,
        vector_length: int,
        matrix_value_semantics: str,
        need_ibs: bool = False,
        gene_ranges=None,
    ) -> tuple[np.ndarray, np.ndarray, None, None, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
            None,
            None,
        )

    def download_stub(compute_backend, totals_tensors, shared_tensors):
        max_len = max((len(tensor) for tensor in totals_tensors), default=0)
        out = np.zeros((len(totals_tensors), 2, max_len), dtype=np.int64)
        for idx, (totals_tensor, shared_tensor) in enumerate(zip(totals_tensors, shared_tensors)):
            out[idx, 0, : len(totals_tensor)] = totals_tensor
            out[idx, 1, : len(shared_tensor)] = shared_tensor
        return out

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "_load_anchor_queue_batch_for_hdf5_torch", tracking_anchor_load)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)
    monkeypatch.setattr(mp, "_prepare_torch_matrix", identity_prepare)
    monkeypatch.setattr(mp, "_compare_anchor_against_target_chunk_torch_device", zero_compare)
    monkeypatch.setattr(mp, "_download_torch_result_tensor_batch", download_stub)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        anchor_queue_size=2,
        backend="torch",
        calculate="ani",
    )

    assert summary.target_chunks == 3
    assert 2 in batch_sizes


def test_matrix_compare_torch_target_queue_prefetches_blocks(tmp_path, monkeypatch):
    torch_module = pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = torch_module

    original_prefetch = mp._load_target_prefetch_unit_for_hdf5_torch
    prefetched_block_sizes: list[int] = []

    def tracking_prefetch(*args, **kwargs):
        block_rows = kwargs.get("block_rows")
        if block_rows is None:
            block_rows = args[3]
        prefetched_block_sizes.append(len(block_rows))
        return original_prefetch(*args, **kwargs)

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    def identity_prepare(compute_backend, matrix, matrix_value_semantics):
        return matrix

    def zero_compare(
        compute_backend,
        anchor_torch,
        target_torch,
        vector_length: int,
        matrix_value_semantics: str,
        need_ibs: bool = False,
        gene_ranges=None,
    ) -> tuple[np.ndarray, np.ndarray, None, None, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
            None,
            None,
        )

    def download_stub(compute_backend, totals_tensors, shared_tensors):
        max_len = max((len(tensor) for tensor in totals_tensors), default=0)
        out = np.zeros((len(totals_tensors), 2, max_len), dtype=np.int64)
        for idx, (totals_tensor, shared_tensor) in enumerate(zip(totals_tensors, shared_tensors)):
            out[idx, 0, : len(totals_tensor)] = totals_tensor
            out[idx, 1, : len(shared_tensor)] = shared_tensor
        return out

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "_load_target_prefetch_unit_for_hdf5_torch", tracking_prefetch)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)
    monkeypatch.setattr(mp, "_prepare_torch_matrix", identity_prepare)
    monkeypatch.setattr(mp, "_compare_anchor_against_target_chunk_torch_device", zero_compare)
    monkeypatch.setattr(mp, "_download_torch_result_tensor_batch", download_stub)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        target_queue_size=2,
        backend="torch",
        calculate="ani",
    )

    assert summary.target_chunks == 3
    assert prefetched_block_sizes
    assert 2 in prefetched_block_sizes


def test_matrix_compare_torch_resumes_after_interruption(tmp_path, monkeypatch):
    torch_module = pytest.importorskip("torch")
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix.h5"
    output_file = tmp_path / "matrix_compare.duckdb"
    _write_many_profiles_same_genome(profile_dir, sample_count=5)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    class FakeTorchBackend:
        def __init__(self, backend: str):
            self.requested = backend
            self.kind = "torch"
            self.device = "cpu"
            self.torch = torch_module

    def two_target_plan(
        vector_length: int,
        remaining_targets: int,
        dtype_name: str,
        memory_limit_bytes: int,
        backend_kind: str,
    ) -> tuple[int, int]:
        return min(2, remaining_targets), vector_length

    def identity_prepare(compute_backend, matrix, matrix_value_semantics):
        return matrix

    def zero_compare(
        compute_backend,
        anchor_torch,
        target_torch,
        vector_length: int,
        matrix_value_semantics: str,
        need_ibs: bool = False,
        gene_ranges=None,
    ) -> tuple[np.ndarray, np.ndarray, None, None, None]:
        return (
            np.zeros(target_torch.shape[2], dtype=np.int64),
            np.zeros(target_torch.shape[2], dtype=np.int64),
            None,
            None,
            None,
        )

    def download_stub(compute_backend, totals_tensors, shared_tensors):
        max_len = max((len(tensor) for tensor in totals_tensors), default=0)
        out = np.zeros((len(totals_tensors), 2, max_len), dtype=np.int64)
        for idx, (totals_tensor, shared_tensor) in enumerate(zip(totals_tensors, shared_tensors)):
            out[idx, 0, : len(totals_tensor)] = totals_tensor
            out[idx, 1, : len(shared_tensor)] = shared_tensor
        return out

    original_mark = mp._mark_completed_pair_genomes
    mark_calls = {"count": 0}

    def interrupt_after_first_commit(compare_conn, completed_pair_genomes):
        mark_calls["count"] += 1
        original_mark(compare_conn, completed_pair_genomes)
        if mark_calls["count"] == 2:
            raise KeyboardInterrupt("simulated ctrl-c")

    monkeypatch.setattr(mp, "MatrixPairComputeBackend", FakeTorchBackend)
    monkeypatch.setattr(mp, "MATRIX_COMPARE_TORCH_CHECKPOINT_BATCH_UNITS", 1)
    monkeypatch.setattr(mp, "_plan_chunk_sizes", two_target_plan)
    monkeypatch.setattr(mp, "_prepare_torch_matrix", identity_prepare)
    monkeypatch.setattr(mp, "_compare_anchor_against_target_chunk_torch_device", zero_compare)
    monkeypatch.setattr(mp, "_download_torch_result_tensor_batch", download_stub)
    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", interrupt_after_first_commit)

    with pytest.raises(KeyboardInterrupt, match="simulated ctrl-c"):
        mp.matrix_compare(
            matrix_db_file=matrix_db,
            output_file=output_file,
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            memory_limit_gb=1.0,
            backend="torch",
            calculate="ani",
        )

    _metadata, completed_pairs_partial, _partial_results = _load_matrix_compare_db(output_file)
    assert completed_pairs_partial == [(0, 1)]

    monkeypatch.setattr(mp, "_mark_completed_pair_genomes", original_mark)
    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="torch",
        calculate="ani",
    )

    assert summary.requested_pairs == 9
    _metadata, completed_pairs, _results = _load_matrix_compare_db(output_file)
    assert completed_pairs == [
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
    ]


def test_build_matrix_db_inserts_separator_rows_for_multiscaffold_genome(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_multi.h5"
    _write_profiles_multiscaffold_same_genome(profile_dir)

    summary = mp.build_matrix_db(
        profile_dir=profile_dir,
        output_file=matrix_db,
        memory_limit_gb=1.0,
    )

    assert summary.sample_count == 2
    _metadata, samples, genomes, scaffolds, matrices, _genes = _load_matrix_hdf5_store(matrix_db)
    assert [sample_name for _sample_idx, sample_name in samples] == ["sample_a", "sample_b"]
    assert [
        (spec.genome_idx, spec.genome, spec.matrix_length, spec.true_length, spec.scaffold_count)
        for spec in genomes
    ] == [(0, "genome1", 5, 4, 2)]
    assert [
        (spec.genome_idx, spec.scaffold_ordinal, spec.genome, spec.chrom, spec.axis_start, spec.axis_end, spec.vector_length)
        for spec in scaffolds
    ] == [
        (0, 0, "genome1", "chr1", 0, 1, 2),
        (0, 1, "genome1", "chr2", 3, 4, 2),
    ]
    # row 2 is the synthetic separator row between the two scaffolds
    assert matrices["0"][0][2].tolist() == [0, 0, 0, 0]
    assert matrices["0"][1][2].tolist() == [0, 0, 0, 0]


def test_matrix_compare_ibs_resets_at_separator_rows(tmp_path):
    profile_dir = tmp_path / "profiles"
    matrix_db = tmp_path / "matrix_multi.h5"
    output_file = tmp_path / "matrix_multi_compare.duckdb"
    _write_profiles_multiscaffold_same_genome(profile_dir)
    mp.build_matrix_db(profile_dir=profile_dir, output_file=matrix_db, memory_limit_gb=1.0)

    summary = mp.matrix_compare(
        matrix_db_file=matrix_db,
        output_file=output_file,
        min_cov=mp.MATRIX_BUILD_MIN_COV,
        memory_limit_gb=1.0,
        backend="numpy",
        calculate="ani+ibs",
    )

    assert summary.requested_pairs == 1
    _metadata, completed_pairs, results = _load_matrix_compare_db(output_file)
    assert completed_pairs == [(0, 1)]
    actual = results.select(
        ["sample_1", "sample_2", "genome", "total_positions", "share_allele_pos", "genome_pop_ani", "max_consecutive_length"]
    ).sort(["sample_1", "sample_2", "genome"])
    expected = (
        cp.compare_genomes(
            mpile_contig_1=profile_dir / "sample_a.parquet",
            mpile_contig_2=profile_dir / "sample_b.parquet",
            min_cov=mp.MATRIX_BUILD_MIN_COV,
            genome_scope="all",
            ani_method="popani",
            engine="polars",
            calculate="ani+ibs",
            stb_file=None,
        )
        .collect(engine="streaming")
        .with_columns(
            sample_1=pl.lit("sample_a"),
            sample_2=pl.lit("sample_b"),
        )
        .select(
            [
                "sample_1",
                "sample_2",
                "genome",
                "total_positions",
                "share_allele_pos",
                "genome_pop_ani",
                "max_consecutive_length",
            ]
        )
        .sort(["sample_1", "sample_2", "genome"])
    )

    assert actual.equals(expected)
    assert actual.get_column("max_consecutive_length").to_list() == [2]
