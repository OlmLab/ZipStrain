import asyncio
import pathlib

import polars as pl
import pytest

from zipstrain import database, task_manager


class WriteOutputTask(task_manager.Task):
    TEMPLATE_CMD = "echo done > <output-file>"


class FailTask(task_manager.Task):
    TEMPLATE_CMD = "echo fail >&2; exit 1"


class SleepTask(task_manager.Task):
    TEMPLATE_CMD = "sleep 3; echo done > <output-file>"


class BrokenPlaceholderTask(task_manager.Task):
    TEMPLATE_CMD = "echo <missing>"


class _SingleTaskGenerator(task_manager.TaskGenerator):
    def __init__(self, task: task_manager.Task):
        self._task = task
        super().__init__(data=None, yield_size=1)

    def get_total_tasks(self) -> int:
        return 1

    async def generate_tasks(self):
        yield [self._task]


class _FinalBatchRunner(task_manager.Runner):
    def __init__(self, run_dir: pathlib.Path):
        final_task = WriteOutputTask(
            id="final_task",
            inputs={},
            expected_outputs={"output-file": task_manager.FileOutput("final_ok.txt")},
            engine=task_manager.LocalEngine(""),
        )
        super().__init__(
            run_dir=run_dir,
            task_generator=_SingleTaskGenerator(
                WriteOutputTask(
                    id="base_task",
                    inputs={},
                    expected_outputs={"output-file": task_manager.FileOutput("base_ok.txt")},
                    engine=task_manager.LocalEngine(""),
                )
            ),
            container_engine=task_manager.LocalEngine(""),
            batch_factory=task_manager.LocalBatch,
            final_batch_factory=task_manager.LocalBatch,
            max_concurrent_batches=1,
            poll_interval=0.05,
            tasks_per_batch=1,
            batch_type="local",
            slurm_config=None,
        )
        self._final_task = final_task

    async def _batcher(self):
        task = await self.tasks_queue.get()
        if task is not None:
            await self.batches_queue.put(
                task_manager.LocalBatch(
                    tasks=[task],
                    id="batch_0",
                    run_dir=self.run_dir,
                    expected_outputs=[],
                )
            )
        await self.batches_queue.put(None)
        self._batcher_done = True

    def _create_final_batch(self):
        if self._final_batch_created:
            return None
        return task_manager.LocalBatch(
            tasks=[self._final_task],
            id="Outputs",
            run_dir=self.run_dir,
            expected_outputs=[],
        )


def test_create_new_input_type():
    class NewInput(task_manager.Input):
        def validate(self) -> None:
            if not isinstance(self.value, str):
                raise ValueError(f"Input value {self.value} is not a string.")

        def get_value(self) -> str:
            return self.value

    with pytest.raises(ValueError, match="Input value 123 is not a string."):
        NewInput(123)


def test_file_input(tmp_path):
    with pytest.raises(FileNotFoundError, match="Input file .* does not exist."):
        task_manager.FileInput("non_existent_file.txt")

    test_file = tmp_path / "test.txt"
    test_file.write_text("Just a test file.")
    file_input = task_manager.FileInput(str(test_file))
    assert file_input.get_value() == str(test_file.absolute())


def test_string_input():
    string_input = task_manager.StringInput("test_string")
    assert string_input.get_value() == "test_string"

    with pytest.raises(ValueError, match="Input value 123 is not a string."):
        task_manager.StringInput(123)


def test_int_input():
    int_input = task_manager.IntInput(42)
    assert int_input.get_value() == "42"

    with pytest.raises(ValueError, match="Input value 'test_string' is not an integer."):
        task_manager.IntInput("test_string")


def test_engine_wrap(tmp_path):
    test_file1 = tmp_path / "file1.txt"
    test_file1.write_text("File 1")
    test_file2 = tmp_path / "file2.txt"
    test_file2.write_text("File 2")

    file_inputs = [task_manager.FileInput(str(test_file1)), task_manager.FileInput(str(test_file2))]
    docker_engine = task_manager.DockerEngine("my_docker_image")
    apptainer_engine = task_manager.ApptainerEngine("my_apptainer_image")
    local_engine = task_manager.LocalEngine("my_local_image")
    command = "echo Hello World"
    wrapped_command = docker_engine.wrap(command, file_inputs)

    expected_mounts = f"-v {test_file1.absolute()}:{test_file1.absolute()} -v {test_file2.absolute()}:{test_file2.absolute()}"
    expected_command = f"docker run {expected_mounts} my_docker_image {command}"

    assert wrapped_command == expected_command
    wrapped_command_apptainer = apptainer_engine.wrap(command, file_inputs)
    expected_mounts_apptainer = f"--bind {test_file1.absolute()}:{test_file1.absolute()},{test_file2.absolute()}:{test_file2.absolute()}"
    expected_command_apptainer = f"apptainer run {expected_mounts_apptainer} my_apptainer_image {command}"
    assert wrapped_command_apptainer == expected_command_apptainer
    wrapped_command_local = local_engine.wrap(command, file_inputs)
    expected_command_local = f"{command}"
    assert wrapped_command_local == expected_command_local


def test_task_map_io_unmapped_placeholder_raises():
    task = BrokenPlaceholderTask(
        id="broken",
        inputs={},
        expected_outputs={},
        engine=task_manager.LocalEngine(""),
    )
    with pytest.raises(ValueError, match="Remaining placeholders"):
        task.map_io()


def test_local_batch_success_writes_logs(tmp_path):
    run_dir = tmp_path / "run"
    task = WriteOutputTask(
        id="task_1",
        inputs={},
        expected_outputs={"output-file": task_manager.FileOutput("ok.txt")},
        engine=task_manager.LocalEngine(""),
    )
    batch = task_manager.LocalBatch(tasks=[task], id="batch_0", run_dir=run_dir, expected_outputs=[])
    batch._set_file_semaphore(asyncio.Semaphore(2))

    asyncio.run(batch.run())

    assert batch.status == task_manager.Status.SUCCESS.value
    assert (run_dir / "batch_0" / "task_1" / "ok.txt").exists()
    assert (run_dir / "batch_0" / ".status").read_text().strip() == task_manager.Status.DONE.value

    run_log = (run_dir / "batch_events.log").read_text()
    batch_log = (run_dir / "batch_0" / "batch.log").read_text()
    assert "BATCH batch_0 START" in run_log
    assert "BATCH batch_0 DONE" in run_log
    assert "progress=1/1" in run_log
    assert "BATCH batch_0" in batch_log


def test_local_batch_runtime_failure_sets_failed_and_logs(tmp_path):
    run_dir = tmp_path / "run"
    task = FailTask(
        id="task_fail",
        inputs={},
        expected_outputs={},
        engine=task_manager.LocalEngine(""),
    )
    batch = task_manager.LocalBatch(tasks=[task], id="batch_fail", run_dir=run_dir, expected_outputs=[])
    batch._set_file_semaphore(asyncio.Semaphore(2))

    with pytest.raises(RuntimeError, match="hit the following error"):
        asyncio.run(batch.run())

    assert batch.status == task_manager.Status.FAILED.value
    assert (run_dir / "batch_fail" / ".status").read_text().strip() == task_manager.Status.FAILED.value
    run_log = (run_dir / "batch_events.log").read_text()
    assert "BATCH batch_fail FAILED" in run_log
    assert "runtime error" in run_log


def test_local_batch_missing_batch_outputs_sets_failed_and_logs(tmp_path):
    run_dir = tmp_path / "run"
    task = WriteOutputTask(
        id="task_ok",
        inputs={},
        expected_outputs={"output-file": task_manager.FileOutput("ok.txt")},
        engine=task_manager.LocalEngine(""),
    )
    batch = task_manager.LocalBatch(
        tasks=[task],
        id="batch_missing_output",
        run_dir=run_dir,
        expected_outputs=[task_manager.BatchFileOutput("required_batch_output.txt")],
    )
    batch._set_file_semaphore(asyncio.Semaphore(2))

    with pytest.raises(FileNotFoundError, match="expected output is missing"):
        asyncio.run(batch.run())

    assert batch.status == task_manager.Status.FAILED.value
    run_log = (run_dir / "batch_events.log").read_text()
    assert "BATCH batch_missing_output FAILED" in run_log
    assert "missing expected outputs" in run_log


def test_local_batch_cancel_sets_failed_and_logs(tmp_path):
    run_dir = tmp_path / "run"
    task = SleepTask(
        id="task_sleep",
        inputs={},
        expected_outputs={"output-file": task_manager.FileOutput("ok.txt")},
        engine=task_manager.LocalEngine(""),
    )
    batch = task_manager.LocalBatch(tasks=[task], id="batch_cancel", run_dir=run_dir, expected_outputs=[])
    batch._set_file_semaphore(asyncio.Semaphore(2))

    async def _run_and_cancel() -> None:
        task_future = asyncio.create_task(batch.run())
        await asyncio.sleep(0.3)
        await batch.cancel()
        with pytest.raises(RuntimeError):
            await task_future

    asyncio.run(_run_and_cancel())
    assert batch.status == task_manager.Status.FAILED.value
    run_log = (run_dir / "batch_events.log").read_text()
    assert "BATCH batch_cancel FAILED" in run_log
    assert task_manager.Messages.CANCELLED_BY_USER.value in run_log


def test_batch_initial_status_reads_success_marker(tmp_path):
    run_dir = tmp_path / "run"
    batch_dir = run_dir / "batch_existing"
    batch_dir.mkdir(parents=True)
    (batch_dir / ".status").write_text(task_manager.Status.SUCCESS.value)
    batch = task_manager.LocalBatch(tasks=[], id="batch_existing", run_dir=run_dir, expected_outputs=[])
    assert batch.status == task_manager.Status.SUCCESS.value


def test_profile_task_generator_includes_gene_stats_output(tmp_path):
    bam = tmp_path / "sample_1.bam"
    bam.write_text("dummy")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    null_model_file = tmp_path / "null_model.parquet"
    null_model_file.write_text("dummy")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t10\n")
    gene_range_file = tmp_path / "gene_info.parquet"
    gene_range_file.write_text("chr1\t0\t10\tgene1\n")
    genome_length_file = tmp_path / "genome_lengths.parquet"
    genome_length_file.write_text("dummy")

    data = pl.DataFrame({"sample_name": ["sample_1"], "bamfile": [str(bam)]}).lazy()
    generator = task_manager.ProfileTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        reference_fasta_file=str(reference_fasta),
        stb_file=str(stb_file),
        null_model_file=str(null_model_file),
        profile_bed_file=str(bed_file),
        gene_range_file=str(gene_range_file),
        profiling_contract_file=None,
        genome_length_file=str(genome_length_file),
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    expected_outputs = tasks[0].expected_outputs
    assert set(expected_outputs.keys()) == {"profile", "genome-stats", "gene-stats"}
    assert expected_outputs["profile"]._expected_file_name == "sample_1_profile.parquet"
    assert expected_outputs["genome-stats"]._expected_file_name == "sample_1_genome_stats.parquet"
    assert expected_outputs["gene-stats"]._expected_file_name == "sample_1_gene_stats.parquet"
    assert tasks[0].inputs["min-freq"].get_value() == "0.01"


def test_profile_task_generator_adds_codon_output_when_dnds_is_prepared(tmp_path):
    files = {
        "bam": tmp_path / "sample.bam",
        "reference": tmp_path / "reference.fna",
        "stb": tmp_path / "reference.stb",
        "null": tmp_path / "null.parquet",
        "bed": tmp_path / "reference.bed",
        "genes": tmp_path / "gene_info.parquet",
        "lengths": tmp_path / "lengths.parquet",
    }
    for path in files.values():
        path.write_text("dummy")
    data = pl.DataFrame(
        {"sample_name": ["sample"], "bamfile": [str(files["bam"])]}
    ).lazy()
    generator = task_manager.ProfileTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        reference_fasta_file=str(files["reference"]),
        stb_file=str(files["stb"]),
        null_model_file=str(files["null"]),
        profile_bed_file=str(files["bed"]),
        gene_range_file=str(files["genes"]),
        profiling_contract_file=None,
        genome_length_file=str(files["lengths"]),
        prepare_dnds=True,
        dnds_memory_limit="512MB",
    )

    async def _collect():
        return [task async for tasks in generator.generate_tasks() for task in tasks]

    task = asyncio.run(_collect())[0]
    assert task.inputs["prepare-dnds-arg"].get_value() == "--prepare-dnds"
    assert task.inputs["dnds-memory-limit"].get_value() == "512MB"
    assert task.expected_outputs["codon-profile"]._expected_file_name == "sample_codon_profile.parquet"


def test_profile_bam_task_template_moves_gene_stats():
    cmd = task_manager.ProfileBamTask.TEMPLATE_CMD
    assert "<reference-fasta-link-cmd>" in cmd
    assert "<reference-fasta-arg>" in cmd
    assert "<profiling-contract-link-cmd>" in cmd
    assert "<profiling-contract-arg>" in cmd
    assert "--min-mapq <min-mapq>" in cmd
    assert "--min-baseq <min-baseq>" in cmd
    assert "--min-freq <min-freq>" in cmd
    assert "<min-read-ani-arg>" in cmd
    assert "--read-inclusion <read-inclusion>" in cmd
    assert "--null-model null_model.parquet" in cmd
    assert "mv input_profile.parquet <sample-name>_profile.parquet" in cmd
    assert "mv input_gene_stats.parquet <sample-name>_gene_stats.parquet" in cmd


def test_profile_task_generator_handles_missing_gene_range(tmp_path):
    bam = tmp_path / "sample_1.bam"
    bam.write_text("dummy")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    null_model_file = tmp_path / "null_model.parquet"
    null_model_file.write_text("dummy")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t10\n")
    genome_length_file = tmp_path / "genome_lengths.parquet"
    genome_length_file.write_text("dummy")

    data = pl.DataFrame({"sample_name": ["sample_1"], "bamfile": [str(bam)]}).lazy()
    generator = task_manager.ProfileTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        reference_fasta_file=str(reference_fasta),
        stb_file=str(stb_file),
        null_model_file=str(null_model_file),
        profile_bed_file=str(bed_file),
        gene_range_file=None,
        profiling_contract_file=None,
        genome_length_file=str(genome_length_file),
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["reference-fasta-link-cmd"].get_value() == f"ln -s {reference_fasta.absolute()} reference.fasta"
    assert tasks[0].inputs["reference-fasta-arg"].get_value() == "--reference-fasta reference.fasta"
    assert tasks[0].inputs["gene-info-table-link-cmd"].get_value() == ""
    assert tasks[0].inputs["gene-info-table-arg"].get_value() == ""
    assert tasks[0].inputs["profiling-contract-link-cmd"].get_value() == ""
    assert tasks[0].inputs["profiling-contract-arg"].get_value() == ""


def test_profile_task_generator_handles_missing_reference_fasta(tmp_path):
    bam = tmp_path / "sample_1.bam"
    bam.write_text("dummy")
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    null_model_file = tmp_path / "null_model.parquet"
    null_model_file.write_text("dummy")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t10\n")
    genome_length_file = tmp_path / "genome_lengths.parquet"
    genome_length_file.write_text("dummy")

    data = pl.DataFrame({"sample_name": ["sample_1"], "bamfile": [str(bam)]}).lazy()
    generator = task_manager.ProfileTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        reference_fasta_file=None,
        stb_file=str(stb_file),
        null_model_file=str(null_model_file),
        profile_bed_file=str(bed_file),
        gene_range_file=None,
        profiling_contract_file=None,
        genome_length_file=str(genome_length_file),
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["reference-fasta-link-cmd"].get_value() == ""
    assert tasks[0].inputs["reference-fasta-arg"].get_value() == ""


def test_profile_task_generator_wires_optional_profiling_contract(tmp_path):
    bam = tmp_path / "sample_1.bam"
    bam.write_text("dummy")
    reference_fasta = tmp_path / "reference.fna"
    reference_fasta.write_text(">chr1\nACGT\n")
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    null_model_file = tmp_path / "null_model.parquet"
    null_model_file.write_text("dummy")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t10\n")
    genome_length_file = tmp_path / "genome_lengths.parquet"
    genome_length_file.write_text("dummy")
    profiling_contract_file = tmp_path / "profiling_contract.json"
    profiling_contract_file.write_text("{}\n")

    data = pl.DataFrame({"sample_name": ["sample_1"], "bamfile": [str(bam)]}).lazy()
    generator = task_manager.ProfileTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        reference_fasta_file=str(reference_fasta),
        stb_file=str(stb_file),
        null_model_file=str(null_model_file),
        profile_bed_file=str(bed_file),
        gene_range_file=None,
        profiling_contract_file=str(profiling_contract_file),
        genome_length_file=str(genome_length_file),
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["profiling-contract-link-cmd"].get_value() == (
        f"ln -s {profiling_contract_file.absolute()} profiling_contract.json"
    )
    assert tasks[0].inputs["profiling-contract-arg"].get_value() == (
        "--profiling-contract profiling_contract.json"
    )


def test_compare_task_generator_creates_tasks_from_profile_locations(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    profile_3 = tmp_path / "sample_3.parquet"
    for path in (profile_1, profile_2, profile_3):
        path.write_text("dummy")

    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1", "s1"],
            "sample_name_2": ["s2", "s3"],
            "profile_location_1": [str(profile_1), str(profile_1)],
            "profile_location_2": [str(profile_2), str(profile_3)],
        }
    ).lazy()
    config = database.GenomeComparisonConfig(
        scope="all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=str(stb_file),
    )
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 2
    assert tasks[0].inputs["profile_1_file"].get_value() == str(profile_1.absolute())
    assert tasks[0].inputs["profile_2_file"].get_value() == str(profile_2.absolute())
    assert tasks[0].inputs["stb-file-option"].get_value() == "--stb-file"
    assert tasks[0].inputs["stb-file"].get_value() == str(stb_file.absolute())
    assert tasks[0].inputs["ani-method-arg"].get_value() == "--ani-method popani"
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate all"
    assert tasks[0].inputs["duckdb-memory-limit-arg"].get_value() == ""
    assert tasks[0].inputs["duckdb-threads-arg"].get_value() == ""
    assert tasks[0].inputs["compare-engine-arg"].get_value() == "--engine polars"
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate all"
    assert tasks[0].expected_outputs["output-file"]._expected_file_name == "s1_s2_comparison.parquet"


def test_compare_task_generator_adds_duckdb_memory_and_threads_args(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    profile_1.write_text("dummy")
    profile_2.write_text("dummy")

    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    config = database.GenomeComparisonConfig(
        scope="all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=str(stb_file),
    )
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
        ani_method="conani, popani,conani",
        calculate="ani",
        duckdb_memory_limit="2GB",
        duckdb_threads=6,
        compare_engine="duckdb",
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["ani-method-arg"].get_value() == "--ani-method conani,popani"
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate ani"
    assert tasks[0].inputs["duckdb-memory-limit-arg"].get_value() == "--duckdb-memory-limit 2GB"
    assert tasks[0].inputs["duckdb-threads-arg"].get_value() == "--duckdb-threads 6"
    assert tasks[0].inputs["compare-engine-arg"].get_value() == "--engine duckdb"
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate ani"


def test_compare_task_generator_accepts_calculate_combo_and_forwards_it(tmp_path):
    # Regression: the batched standard-compare path used to reject any
    # --calculate value other than {"all", "ani"}, even though the per-pair
    # single-compare worker understands '+'/','-joined combos.
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    profile_1.write_text("dummy")
    profile_2.write_text("dummy")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    config = database.GenomeComparisonConfig(
        scope="all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=None,
    )
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
        calculate="genome_ani+ibs",
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate genome_ani+ibs"

    with pytest.raises(ValueError):
        task_manager.CompareTaskGenerator(
            data=data,
            yield_size=1,
            container_engine=task_manager.LocalEngine(""),
            comp_config=config,
            calculate="ani+bogus",
        )


def test_compare_task_generator_omits_stb_arg_when_not_configured(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    profile_1.write_text("dummy")
    profile_2.write_text("dummy")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    config = database.GenomeComparisonConfig(
        scope="all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=None,
    )
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["stb-file-option"].get_value() == ""
    assert tasks[0].inputs["stb-file"].get_value() == ""


def test_fast_compare_batch_cleanup_is_idempotent(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    for path in (profile_1, profile_2):
        path.write_text("dummy")

    task = task_manager.FastCompareTask(
        id="cmp",
        inputs={
            "profile_1_file": task_manager.FileInput(profile_1),
            "profile_2_file": task_manager.FileInput(profile_2),
            "stb-file-option": task_manager.StringInput(""),
            "stb-file": task_manager.StringInput(""),
            "gene-info-table-option": task_manager.StringInput(""),
            "gene-info-table": task_manager.StringInput(""),
            "min_cov": task_manager.IntInput(5),
            "min-gene-compare-len": task_manager.IntInput(100),
            "ani-method-arg": task_manager.StringInput("--ani-method popani"),
            "calculate-arg": task_manager.StringInput("--calculate ani"),
            "duckdb-memory-limit-arg": task_manager.StringInput(""),
            "duckdb-threads-arg": task_manager.StringInput(""),
            "compare-engine-arg": task_manager.StringInput("--engine polars"),
            "scope": task_manager.StringInput("all"),
        },
        expected_outputs={"output-file": task_manager.FileOutput("out.parquet")},
        engine=task_manager.LocalEngine(""),
    )
    batch = task_manager.FastCompareLocalBatch(tasks=[task], id="batch_0", run_dir=tmp_path, expected_outputs=[])
    batch.cleanup()
    batch.cleanup()
    assert batch._cleaned_up is True


def test_unified_compare_task_generator_adds_gene_and_duckdb_args(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    profile_1.write_text("dummy")
    profile_2.write_text("dummy")

    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    gene_range_table = tmp_path / "gene_info.parquet"
    gene_range_table.write_text("gene1\tchr1\t1\t10\n")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    config = database.GenomeComparisonConfig(
        gene_range_table_loc=str(gene_range_table),
        scope="all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=str(stb_file),
    )
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
        duckdb_memory_limit="3GB",
        duckdb_threads=4,
        compare_engine="duckdb",
        calculate="gene",
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["stb-file-option"].get_value() == "--stb-file"
    assert tasks[0].inputs["stb-file"].get_value() == str(stb_file.absolute())
    assert tasks[0].inputs["duckdb-memory-limit-arg"].get_value() == "--duckdb-memory-limit 3GB"
    assert tasks[0].inputs["duckdb-threads-arg"].get_value() == "--duckdb-threads 4"
    assert tasks[0].inputs["compare-engine-arg"].get_value() == "--engine duckdb"
    assert tasks[0].inputs["gene-info-table-option"].get_value() == "--gene-info-table"
    assert tasks[0].inputs["gene-info-table"].get_value() == str(gene_range_table.absolute())
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate gene"


def test_compare_task_generator_stages_codon_sidecars_for_dnds(tmp_path):
    profile_1 = tmp_path / "sample_1_profile.parquet"
    profile_2 = tmp_path / "sample_2_profile.parquet"
    codon_1 = tmp_path / "sample_1_codon_profile.parquet"
    codon_2 = tmp_path / "sample_2_codon_profile.parquet"
    gene_info = tmp_path / "gene_info.parquet"
    for path in (profile_1, profile_2, codon_1, codon_2, gene_info):
        path.write_text("dummy")
    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=database.GenomeComparisonConfig(
            gene_range_table_loc=str(gene_info),
            scope="all",
            min_cov=5,
            min_gene_compare_len=1,
            stb_file_loc=None,
        ),
        calculate="gene",
        dnds=True,
        dnds_min_major_freq=0.8,
    )

    async def _collect():
        return [task async for tasks in generator.generate_tasks() for task in tasks]

    task = asyncio.run(_collect())[0]
    assert task.inputs["dnds-arg"].get_value() == "--dnds --dnds-min-major-freq 0.8"
    assert task.inputs["codon-profile-1-option"].get_value() == "--codon-profile-1"
    assert task.inputs["codon-profile-1"].get_value() == str(codon_1.absolute())
    assert task.inputs["codon-profile-2"].get_value() == str(codon_2.absolute())


def test_unified_gene_task_generator_omits_stb_arg_when_not_configured(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    profile_1.write_text("dummy")
    profile_2.write_text("dummy")
    gene_range_table = tmp_path / "gene_info.parquet"
    gene_range_table.write_text("gene1\tchr1\t1\t10\n")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    config = database.GenomeComparisonConfig(
        gene_range_table_loc=str(gene_range_table),
        scope="all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=None,
    )
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
        calculate="gene",
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["stb-file-option"].get_value() == ""
    assert tasks[0].inputs["stb-file"].get_value() == ""


def test_unified_compare_container_mounts_stb_and_gene_ranges(tmp_path):
    profile_1 = tmp_path / "sample_1.parquet"
    profile_2 = tmp_path / "sample_2.parquet"
    stb_file = tmp_path / "stb.tsv"
    gene_range_table = tmp_path / "gene_info.parquet"
    for path in (profile_1, profile_2, stb_file, gene_range_table):
        path.write_text("dummy")

    data = pl.DataFrame(
        {
            "sample_name_1": ["s1"],
            "sample_name_2": ["s2"],
            "profile_location_1": [str(profile_1)],
            "profile_location_2": [str(profile_2)],
        }
    ).lazy()
    engine = task_manager.DockerEngine("zipstrain:test")
    generator = task_manager.CompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=engine,
        comp_config=database.GenomeComparisonConfig(
            scope="all",
            min_cov=5,
            min_gene_compare_len=100,
            stb_file_loc=str(stb_file),
            gene_range_table_loc=str(gene_range_table),
        ),
        calculate="gene",
    )

    async def _collect():
        async for chunk in generator.generate_tasks():
            return chunk

    tasks = asyncio.run(_collect())
    task_manager.FastCompareLocalBatch(
        tasks=tasks,
        id="batch_0",
        run_dir=tmp_path / "run",
        expected_outputs=[],
    )
    command = tasks[0].command
    for path in (profile_1, profile_2, stb_file, gene_range_table):
        absolute = path.absolute()
        assert f"-v {absolute}:{absolute}" in command


def test_fast_compare_templates_call_utilities_single_compare():
    assert "zipstrain utilities single-compare" in task_manager.FastCompareTask.TEMPLATE_CMD
    assert "<stb-file-option> <stb-file>" in task_manager.FastCompareTask.TEMPLATE_CMD
    assert (
        "<gene-info-table-option> <gene-info-table>"
        in task_manager.FastCompareTask.TEMPLATE_CMD
    )


def test_slurm_config_validation_and_args():
    conf = task_manager.SlurmConfig(
        time="10:00:00",
        tasks=2,
        mem=16,
        additional_params={"cpus-per-task": "4"},
    )
    args = conf.to_slurm_args()
    assert "#SBATCH --time=10:00:00" in args
    assert "#SBATCH --ntasks=2" in args
    assert "#SBATCH --mem=16G" in args
    assert "#SBATCH --cpus-per-task=4" in args

    with pytest.raises(ValueError, match="Time must be in the format"):
        task_manager.SlurmConfig(time="1:2", tasks=1, mem=1)


def test_runner_exits_after_final_batch(tmp_path):
    run_dir = tmp_path / "run"
    runner = _FinalBatchRunner(run_dir=run_dir)
    asyncio.run(asyncio.wait_for(runner.run(), timeout=10))
    assert (run_dir / "batch_0" / "base_task" / "base_ok.txt").exists()
    assert (run_dir / "Outputs" / "final_task" / "final_ok.txt").exists()


def test_collect_comps_template_is_retry_safe():
    cmd = task_manager.CollectComps.TEMPLATE_CMD
    assert "rm -rf comps" in cmd
    assert '! -path "./comps/*"' in cmd
    assert "cp */*_comparison.parquet comps/" not in cmd


def _make_fake_batch_run(run_dir: pathlib.Path, sample_names, batches):
    """Build a fake raw profile-run tree: run_dir/batch_N/<sample>/ with outputs
    plus intermediate files, batch-level logs, and a top-level batch_events.log."""
    (run_dir / task_manager.Batch.RUN_LOG_FILE).write_text("run log\n")
    for batch_idx, batch_samples in enumerate(batches):
        batch_dir = run_dir / f"batch_{batch_idx}"
        batch_dir.mkdir(parents=True)
        (batch_dir / f"batch_{batch_idx}.sh").write_text("#!/bin/bash\n")
        (batch_dir / f"batch_{batch_idx}.out").write_text("stdout\n")
        (batch_dir / f"batch_{batch_idx}.err").write_text("")
        (batch_dir / ".status").write_text("done\n")
        (batch_dir / "batch.log").write_text("batch log\n")
        for sample_name in batch_samples:
            sample_dir = batch_dir / sample_name
            sample_dir.mkdir()
            # Real outputs
            for suffix in ("_profile.parquet", "_genome_stats.parquet", "_gene_stats.parquet"):
                (sample_dir / f"{sample_name}{suffix}").write_text("data")
            # Intermediates
            (sample_dir / ".status").write_text("done\n")
            (sample_dir / "input.bam").write_text("bam")
            (sample_dir / "bed_file.bed").write_text("bed")
    return sample_names


def test_reorganize_profile_run_output_flattens_and_tidies(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    samples = _make_fake_batch_run(run_dir, ["s1", "s2", "s3"], batches=[["s1", "s2"], ["s3"]])

    task_manager._reorganize_profile_run_output(run_dir)

    # No batch directories remain.
    assert list(run_dir.glob("batch_*")) == []

    # Batch logs collected under profiling_assets/log.
    log_dir = run_dir / "profiling_assets" / "log"
    assert (log_dir / "batch_events.log").exists()
    assert (log_dir / "batch_0.out").exists()
    assert (log_dir / "batch_1.out").exists()
    # Colliding names (per-batch .status / batch.log) are all kept: the first
    # keeps its name, later collisions are prefixed by their batch directory.
    assert (log_dir / ".status").exists()
    assert (log_dir / "batch_1_.status").exists()
    assert (log_dir / "batch.log").exists()
    assert (log_dir / "batch_1_batch.log").exists()

    for sample_name in samples:
        sample_dir = run_dir / sample_name
        top_level = {p.name for p in sample_dir.iterdir()}
        assert top_level == {
            f"{sample_name}_profile.parquet",
            f"{sample_name}_genome_stats.parquet",
            f"{sample_name}_gene_stats.parquet",
            "intermediate_files",
        }
        intermediates = {p.name for p in (sample_dir / "intermediate_files").iterdir()}
        assert intermediates == {".status", "input.bam", "bed_file.bed"}


def _make_profile_runner(run_dir):
    generator = _SingleTaskGenerator(
        WriteOutputTask(
            id="s1",
            inputs={},
            expected_outputs={"output-file": task_manager.FileOutput("s1_ok.txt")},
            engine=task_manager.LocalEngine(""),
        )
    )
    return task_manager.ProfileRunner(
        run_dir=run_dir,
        task_generator=generator,
        container_engine=task_manager.LocalEngine(""),
    )


def test_final_summary_reports_batch_count_elapsed_and_paths(tmp_path):
    run_dir = tmp_path / "run"
    runner = _make_profile_runner(run_dir)
    runner._success_batches_count = 3
    runner._failed_batches_count = 0
    runner._produced_tasks_count = 5

    text = runner._final_summary_text(total_batches=3, elapsed_str="0:01:07")

    assert "Run finished!" in text
    # Correct batch fraction (not the old total_tasks/tasks_per_batch bug).
    assert "3/3 batches succeeded." in text
    assert "Produced tasks: 5" in text
    assert "Elapsed: 0:01:07" in text
    assert f"Output: {run_dir.absolute()}" in text
    assert f"Logs:   {(run_dir / 'profiling_assets' / 'log').absolute()}" in text


def test_final_summary_reports_failures(tmp_path):
    runner = _make_profile_runner(tmp_path / "run")
    runner._success_batches_count = 2
    runner._failed_batches_count = 1
    runner._produced_tasks_count = 3

    text = runner._final_summary_text(total_batches=3, elapsed_str="0:00:30")

    assert "Run finished with failures." in text
    assert "2/3 batches succeeded (1 failed)." in text


def _make_fake_compare_run(run_dir, batch_name, output_filename):
    """Build a fake raw compare-run tree: a work batch + an Outputs final batch."""
    (run_dir / task_manager.Batch.RUN_LOG_FILE).write_text("run log\n")

    work = run_dir / batch_name
    (work / "concat_parquet").mkdir(parents=True)
    (work / f"{batch_name}.err").write_text("")
    (work / f"{batch_name}.out").write_text("out")
    (work / f"{batch_name}.sh").write_text("#!/bin/bash\n")
    (work / "batch.log").write_text("log")
    (work / ".status").write_text("done\n")
    (work / "concat_parquet" / ".status").write_text("done\n")
    (work / "concat_parquet" / f"Merged_{batch_name}.parquet").write_text("data")

    outputs = run_dir / "Outputs"
    (outputs / "prepare_outputs").mkdir(parents=True)
    (outputs / output_filename).write_text("MERGED")  # the real output
    (outputs / "Outputs.err").write_text("")
    (outputs / "Outputs.out").write_text("out")
    (outputs / "Outputs.sh").write_text("#!/bin/bash\n")
    (outputs / "batch.log").write_text("log")
    (outputs / ".status").write_text("done\n")
    (outputs / "prepare_outputs" / ".status").write_text("done\n")


def test_reorganize_compare_run_output_flattens_and_tidies(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    batch_name = "batch_0"
    output_filename = "all_comparisons.parquet"
    _make_fake_compare_run(run_dir, batch_name, output_filename)

    task_manager._reorganize_compare_run_output(run_dir)

    # The real merged output is lifted to the top of run_dir.
    assert (run_dir / output_filename).read_text() == "MERGED"

    # No batch-like directories remain at the top level.
    assert not (run_dir / batch_name).exists()
    assert not (run_dir / "Outputs").exists()

    # Logs collected under run_dir/log (with collisions prefixed).
    log_dir = run_dir / "log"
    assert (log_dir / "batch_events.log").exists()
    assert (log_dir / f"{batch_name}.out").exists()
    assert (log_dir / "Outputs.out").exists()
    assert (log_dir / "batch.log").exists()
    assert (log_dir / "Outputs_batch.log").exists()

    # Intermediate task dirs preserved under intermediate_files/<batch>/.
    inter = run_dir / "intermediate_files"
    assert (inter / batch_name / "concat_parquet" / f"Merged_{batch_name}.parquet").exists()
    assert (inter / "Outputs" / "prepare_outputs" / ".status").exists()


def test_unified_compare_config_carries_gene_range_table(tmp_path):
    """The single task path carries the table needed by gene calculations."""
    gene_range_table = tmp_path / "gene_info.parquet"
    pl.DataFrame(
        {"gene": ["g1"], "scaffold": ["chr1"], "start": [1], "end": [10]}
    ).write_parquet(gene_range_table)

    config = database.GenomeComparisonConfig(
        scope="all",
        min_cov=5,
        min_gene_compare_len=1,
        gene_range_table_loc=str(gene_range_table),
    )
    assert config.gene_range_table_loc == str(gene_range_table)
    assert (
        "<gene-info-table-option> <gene-info-table>"
        in task_manager.FastCompareTask.TEMPLATE_CMD
    )
