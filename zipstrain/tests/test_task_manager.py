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
    stb_file = tmp_path / "stb.tsv"
    stb_file.write_text("chr1\tgenome1\n")
    null_model_file = tmp_path / "null_model.parquet"
    null_model_file.write_text("dummy")
    bed_file = tmp_path / "genomes.bed"
    bed_file.write_text("chr1\t0\t10\n")
    gene_range_file = tmp_path / "gene_ranges.tsv"
    gene_range_file.write_text("chr1\t0\t10\tgene1\n")
    genome_length_file = tmp_path / "genome_lengths.parquet"
    genome_length_file.write_text("dummy")

    data = pl.DataFrame({"sample_name": ["sample_1"], "bamfile": [str(bam)]}).lazy()
    generator = task_manager.ProfileTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        stb_file=str(stb_file),
        null_model_file=str(null_model_file),
        profile_bed_file=str(bed_file),
        gene_range_file=str(gene_range_file),
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
    assert expected_outputs["profile"]._expected_file_name == "sample_1.parquet"
    assert expected_outputs["genome-stats"]._expected_file_name == "sample_1_genome_stats.parquet"
    assert expected_outputs["gene-stats"]._expected_file_name == "sample_1_gene_stats.parquet"


def test_profile_bam_task_template_moves_gene_stats():
    cmd = task_manager.ProfileBamTask.TEMPLATE_CMD
    assert "--null-model null_model.parquet" in cmd
    assert "mv input_gene_stats.parquet <sample-name>_gene_stats.parquet" in cmd


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
        gene_db_id="gene_ref",
        reference_id="ref",
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
    assert tasks[0].inputs["mpile_1_file"].get_value() == str(profile_1.absolute())
    assert tasks[0].inputs["mpile_2_file"].get_value() == str(profile_2.absolute())
    assert tasks[0].inputs["ani-method-arg"].get_value() == "--ani-method popani"
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate all"
    assert tasks[0].inputs["duckdb-memory-limit-arg"].get_value() == ""
    assert tasks[0].inputs["duckdb-threads-arg"].get_value() == ""
    assert tasks[0].inputs["compare-engine-arg"].get_value() == "--engine polars"
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
        gene_db_id="gene_ref",
        reference_id="ref",
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
        ani_method="conani",
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
    assert tasks[0].inputs["ani-method-arg"].get_value() == "--ani-method conani"
    assert tasks[0].inputs["calculate-arg"].get_value() == "--calculate ani"
    assert tasks[0].inputs["duckdb-memory-limit-arg"].get_value() == "--duckdb-memory-limit 2GB"
    assert tasks[0].inputs["duckdb-threads-arg"].get_value() == "--duckdb-threads 6"
    assert tasks[0].inputs["compare-engine-arg"].get_value() == "--engine duckdb"


def test_gene_compare_task_generator_adds_duckdb_memory_and_threads_args(tmp_path):
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
    config = database.GeneComparisonConfig(
        gene_db_id="gene_ref",
        reference_genome_id="ref",
        scope="all:all",
        min_cov=5,
        min_gene_compare_len=100,
        stb_file_loc=str(stb_file),
    )
    generator = task_manager.GeneCompareTaskGenerator(
        data=data,
        yield_size=1,
        container_engine=task_manager.LocalEngine(""),
        comp_config=config,
        duckdb_memory_limit="3GB",
        duckdb_threads=4,
        compare_engine="duckdb",
    )

    async def _collect():
        out = []
        async for chunk in generator.generate_tasks():
            out.extend(chunk)
        return out

    tasks = asyncio.run(_collect())
    assert len(tasks) == 1
    assert tasks[0].inputs["duckdb-memory-limit-arg"].get_value() == "--duckdb-memory-limit 3GB"
    assert tasks[0].inputs["duckdb-threads-arg"].get_value() == "--duckdb-threads 4"
    assert tasks[0].inputs["compare-engine-arg"].get_value() == "--engine duckdb"


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


def test_collect_gene_comps_template_is_retry_safe():
    cmd = task_manager.CollectGeneComps.TEMPLATE_CMD
    assert "rm -rf gene_comps" in cmd
    assert '! -path "./gene_comps/*"' in cmd
    assert "cp */*_gene_comparison.parquet gene_comps/" not in cmd
