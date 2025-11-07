"""zipstrain.task_manager
========================
Lightweight, asyncio-driven orchestration primitives for building and running
scientific data-processing pipelines. This module provides a small, composable
framework for defining Tasks with explicit inputs/outputs, bundling Tasks into
Batches (local or Slurm), and coordinating their execution with a live terminal
UI. It is designed to be easy to extend for new Task types and execution
environments. For most users, this module is not directly used. However, 
it can be used to define new pipelines that chain together multiple steps with clear input/outputs.
The unit of execution is a batch, which is a collection of tasks to be executed together.
Each batch can have an optional finalization step that runs after all tasks are complete.


Key concepts
------------

- Inputs and Outputs:
    These classes encapsulate task inputs and outputs with validation logics.
    By default, Input and output classes for files, strings, and integers are provided.
    If needed, new types can be defined by subclassing Input or Output.
    
- Engines:
    Any task object can use a container engine (Docker or Apptainer) or run natively (LocalEngine).
    
- Task 
    Each task runs a unit of bash script with defined inputs and expected outputs. If an engine is provided,
    the command will be wrapped accordingly to run inside the container.
    
- Batches:
    A batch is a collection of tasks to be executed together. Batches can be run locally or submitted to Slurm.
    Each batch monitors the status of its tasks and updates its own status accordingly. A batch can also have
    expected outputs that are checked after all tasks are complete. Additionally, a batch can have a finalization step that runs after all tasks are complete.
    
- Runner:
    The Runner class orchestrates task generation, batching, and execution. It manages concurrent batch execution,
    monitors progress, and provides a live terminal UI using the rich library.
    
"""


from __future__ import annotations
from enum import StrEnum
import re
import pathlib
from abc import ABC, abstractmethod
import asyncio
import subprocess
import aiofiles
from pydantic import BaseModel, Field, field_validator
from rich.live import Live
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich.console import Console, Group
from rich.panel import Panel
from rich.align import Align
from zipstrain import database
from rich.columns import Columns
import polars as pl
import psutil
import shutil
import signal 


class SlurmConfig(BaseModel):
    """Configuration model for Slurm batch jobs.
    
    Attributes:
        time (str): Time limit for the job in HH:MM:SS format.
        tasks (int): Number of tasks.
        mem (int): Memory in GB.
        additional_params (dict): Additional SLURM parameters as key-value pairs.
        
    NOTE: Additional paramters for slurm should be provided in the additional_params dict in the form
    of {"param-name": "param-value"}, e.g., {"cpus-per-task": "4"} will result in the addition of
    "#SBATCH --cpus-per-task=4" to the sbatch script.
    
    """
    time: str = Field(description="Time limit for the job.")
    tasks: int = Field(default=1, description="Number of tasks.")
    mem: int = Field(default=4, description="Memory in GB.")
    additional_params: dict[str, str] = Field(default_factory=dict, description="Additional SLURM parameters as key-value pairs.")

    @field_validator("time")
    def validate_time(cls, v):
        """Validate time format HH:MM:SS (H..HHH allowed)."""
        if not re.match(r"^\d{1,3}:\d{2}:\d{2}$", v):
            raise ValueError("Time must be in the format HH:MM:SS (H..HHH allowed)")
        return v

    def to_slurm_args(self) -> str:
        """Generates the slurm batch file header form the configuration object"""
        args = [
            f"#SBATCH --time={self.time}",
            f"#SBATCH --ntasks={self.tasks}",
            f"#SBATCH --mem={self.mem}G",
        ]
        for key, value in self.additional_params.items():
            args.append(f"#SBATCH --{key}={value}")
        return "\n".join(args)

    @classmethod
    def from_json(cls, json_path: str | pathlib.Path) -> SlurmConfig:
        """Load SlurmConfig from a JSON file."""
        path = pathlib.Path(json_path)
        if not path.exists():
            raise FileNotFoundError(f"Slurm config file {json_path} does not exist.")
        return cls.model_validate_json(path.read_text())

async def write_file(path: pathlib.Path, text: str, file_semaphore: asyncio.Semaphore) -> None:
    async with file_semaphore:
        async with aiofiles.open(path, "w") as f:
            await f.write(text)

async def read_file(path: pathlib.Path, file_semaphore: asyncio.Semaphore) -> str:
    async with file_semaphore:
        async with aiofiles.open(path, "r") as f:
            content = await f.read()
    return content

class Status(StrEnum):
    """Enumeration of possible task and batch statuses."""
    BATCH_NOT_ASSIGNED = "batch_not_assigned"
    NOT_STARTED = "not_started"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"
    SUBMITTED = "submitted"
    SUCCESS = "success"
    PENDING = "pending"


class Input(ABC):
    """Abstract base class for task inputs. DO NOT INSTANTIATE DIRECTLY.
    Most commonly used Input types are provided but if you want to define a new one,
    subclass this and implement validate() and get_value().
    """
    def __init__(self, value: str | int) -> None:
        self.value = value
        self.validate()

    @abstractmethod
    def validate(self) -> None:
        ...

    @abstractmethod
    def get_value(self) -> str | int:
        ...


class FileInput(Input):
    """This is used when the input is a file path. By default, the validate method checks for file existence."""
    def validate(self, check_exists: bool = True) -> None:
        if check_exists and not pathlib.Path(self.value).exists():
            raise FileNotFoundError(f"Input file {self.value} does not exist.")

    def get_value(self) -> str:
        """Returns the absolute path of the input file as a string."""
        return str(pathlib.Path(self.value).absolute())


class StringInput(Input):
    """This is used when the input is a string."""

    def validate(self) -> None:
        """Validate that the input value is a string."""
        if not isinstance(self.value, str):
            raise ValueError(f"Input value {self.value!r} is not a string.")

    def get_value(self) -> str:
        """Returns the string value."""
        return str(self.value)


class IntInput(Input):
    """This is used when the input is an integer."""
    def validate(self) -> None:
        """Validate that the input value is an integer."""
        if not isinstance(self.value, int):
            raise ValueError(f"Input value {self.value!r} is not an integer.")

    def get_value(self) -> str:
        """
        Returns the integer value as a string.
        """
        return str(self.value)


class Output(ABC):
    """Abstract base class for task outputs. DO NOT INSTANTIATE DIRECTLY.
    Most commonly used Output types are provided but if you want to define a new one,
    subclass this and implement ready().
    This method is used to check if the output is ready/valid after task completion.
    """
    def __init__(self) -> None:
        self._value = None ## Will be set by the task when it completes
        self.task = None  ## Will be set when the output is registered to a task

    @property
    def value(self):
        return self._value

    @abstractmethod
    def ready(self) -> bool:
        ...

    def register_task(self, task: Task) -> None:
        """Registers the task that produces this output. In most cases, you won't need to override this.
        
        Args:
        task (Task): The task that produces this output.
        """
        self.task = task

class FileOutput(Output):
    """This is used when the output is a file path.
    
    Args:
        expected_file (str): The expected output file name relative to the task directory.
    """
    def __init__(self, expected_file:str) -> None:
        self._expected_file_name = expected_file ### When the task is finished, the expected file should be in task.task_dir / expected_file otherwise ready() will return False

    def ready(self) -> bool:
        """Check if the expected output file exists."""
        return True if self.expected_file.absolute().exists() else False

    def register_task(self, task: Task) -> None:
        """Registers the task that produces this output and sets the expected file path.
        
        Args:
        task (Task): The task that produces this output.
        """
        super().register_task(task)
        self.expected_file = self.task.task_dir / self._expected_file_name


class BatchFileOutput(Output):
    """This is used when the output is a file path relative to the batch directory.
    Also it will be registered to the batch instead of the task.
    """
    def __init__(self, expected_file:str) -> None:
        self._expected_file_name = expected_file

    def ready(self) -> bool:
        """Check if the expected output file exists."""
        return True if self.expected_file.absolute().exists() else False

    def register_batch(self, batch: Batch) -> None:
        """Registers the batch that produces this output and sets the expected file path.
        
        Args:
        batch (Batch): The batch that produces this output and sets the expected file path.
        """
        self.expected_file = batch.batch_dir / self._expected_file_name


class StringOutput(Output):
    """This is used when the output is a string."""
    def ready(self) -> bool:
        """Check if the output value is a string."""
        if isinstance(self._value, str):
            return True
        elif self._value is not None:
            raise ValueError(f"Output value for task {self.task.id} is not a string.")
        else:
            return False


class IntOutput(Output):
    """This is used when the output is an integer."""
    def ready(self) -> bool:
        """Check if the output value is an integer."""
        if isinstance(self._value, int):
            return True
        elif self._value is not None:
            raise ValueError(f"Output value for task {self.task.id} is not an integer.")
        else:
            return False
        return False


class Engine(ABC):
    def __init__(self, address: str) -> None:
        self.address = address

    @abstractmethod
    def wrap(self, command: str, file_inputs: list[FileInput]) -> str:
        ...


class DockerEngine(Engine):
    def wrap(self, command: str, file_inputs: list[FileInput]) -> str:
        volume_mounts = " ".join(
            [f"-v {file_input.get_value()}:{file_input.get_value()}" for file_input in file_inputs]
        )
        return f"docker run {volume_mounts} {self.address} {command}"


class ApptainerEngine(Engine):
    def wrap(self, command: str, file_inputs: list[FileInput]) -> str:
        volume_mounts = "--bind " + ",".join(
            [f"{file_input.get_value()}:{file_input.get_value()}" for file_input in file_inputs]
        )
        return f"apptainer run {volume_mounts} {self.address} {command}"


class LocalEngine(Engine):
    def wrap(self, command: str, file_inputs: list[FileInput]) -> str:
        return command


class Task(ABC):
    """Abstract base class for tasks. DO NOT INSTANTIATE DIRECTLY. Any new task type should subclass this
    and implement the TEMPLATE_CMD class attribute. Inputs and expected outputs are specified using <>.
    As an example, if a task has an input file called "input-file" and an expected output file called "output-file",
    the TEMPLATE_CMD could be something like:
        TEMPLATE_CMD = "some_command --input <input-file> --output <output-file>"
    the outputs and inputs will be mapped to the command when map_io() is called later in the runtime.
    """
    TEMPLATE_CMD = ""

    def __init__(
        self,
        id: str,
        inputs: dict[str, Input | Output],
        expected_outputs: dict[str, Output] ,
        engine: Engine,
        batch_obj: Batch | None = None,
        file_semaphore: asyncio.Semaphore | None = None
    ) -> None:
        self.id = id
        self.inputs = inputs
        self.expected_outputs = expected_outputs
        self._batch_obj = batch_obj
        self.engine = engine
        self._status = self._get_initial_status()
        self.file_semaphore = file_semaphore

    def map_io(self) -> None:
        """Maps inputs and expected outputs to the command template. Note that when this method is called,
        all of the inputs and outputs in the TEMPLATE_CMD must be defined in the inputs and expected_outputs dictionaries.
        However, this method is not called by the user directly. It is called by the Batch when the task is added to a batch.
        """
        cmd = self.TEMPLATE_CMD
        for key, value in self.inputs.items():
            cmd = cmd.replace(f"<{key}>", value.get_value())
        # if any placeholders remain, report them

        for handle, output in self.expected_outputs.items():
            cmd = cmd.replace(f"<{handle}>", str(output.expected_file.absolute()))
        remaining = re.findall(r"<\w+>", cmd)
        if remaining:
            raise ValueError(f"Not all inputs were mapped in task {self.id}. Remaining placeholders: {remaining}")
        self._command = cmd

    @property
    def batch_dir(self) -> pathlib.Path:
        """Returns the batch directory path. Raises an error if the task is not associated with any batch yet."""
        if self._batch_obj is None:
            raise ValueError(f"Task {self.id} is not associated with any batch yet.")
        return self._batch_obj.batch_dir

    @property
    def task_dir(self) -> pathlib.Path:
        """Returns the task directory path."""
        return self.batch_dir / self.id
    
    @property
    def command(self) -> str:
        """Returns the command to be executed, wrapped with the engine if applicable."""
        file_inputs = [v for v in self.inputs.values() if isinstance(v, FileInput)]
        return self.engine.wrap(self._command, file_inputs)

    @property
    def pre_run(self) -> str:
        """Does the necessary setup before running the task command. This should not be overridden by subclasses unless a task needs special setup like
        batch aggregation."""
        return f"echo {Status.RUNNING.value} > {self.task_dir.absolute()}/.status && cd {self.task_dir.absolute()}"

    @property
    def status(self) -> str:
        """Returns the current status of the task."""
        return self._status

    @property
    def post_run(self) -> str:
        """Does the necessary steps after running the task command. This should not be overridden by subclasses unless a task needs special teardown like
        batch aggregation."""
        return f"cd {self.batch_dir.absolute()} && echo {Status.DONE.value} > {self.task_dir.absolute()}/.status"

    async def get_status(self) -> str:
        """Asynchronously reads the task status from the .status file in the task directory."""
        status_path = self.task_dir / ".status"
        # read the status file if it exists
        if status_path.exists():
            raw = await read_file(status_path, self.file_semaphore)
            self._status = raw.strip()

            # if task reported 'done', check outputs to decide success/failure
            if self._status == Status.DONE.value:
                all_ready = True
                try:
                    for output in self.expected_outputs.values():
                        if not output.ready():
                            all_ready = False
                            break
                except Exception:
                    all_ready = False

                if all_ready:
                    self._status = Status.SUCCESS.value
                    await write_file(status_path, Status.SUCCESS.value, self.file_semaphore)
                else:
                    self._status = Status.FAILED.value
                    await write_file(status_path, Status.FAILED.value, self.file_semaphore)
                    raise ValueError(f"Task {self.id} reported done but outputs are not ready or invalid. {self.expected_outputs['output-file'].expected_file.absolute()}")

        return self._status

    def _get_initial_status(self) -> str:
        """Returns the initial status of the task based on the presence of the batch and task directories."""
        if self._batch_obj is None:
            return Status.BATCH_NOT_ASSIGNED.value
        if not self.task_dir.exists():
            return Status.NOT_STARTED.value
        status_file = self.task_dir / ".status"
        with open(status_file, mode="r") as f:
            status_as_written = f.read().strip()
        if status_as_written in (Status.DONE.value, Status.SUCCESS.value):
            all_ready = True
            try:
                for output in self.expected_outputs.values():
                    if not output.ready():
                        all_ready = False
                        break
            except Exception:
                all_ready = False

            if all_ready:
                return Status.SUCCESS.value
            else:
                return Status.FAILED.value

class TaskGenerator(ABC):
    """Abstract base class for task generators. DO NOT INSTANTIATE DIRECTLY. A subclass of this class 
    should provide an async generator method called generate_tasks() that yields lists of Task objects in an async manner.
    Some important concepts:
    
    - generate_tasks() is an async generator that yields lists of Task objects.
    
    - yield_size determines how many tasks are generated and yielded at a time.
    
    - get_total_tasks() returns the total number of tasks that can be generated.

    """
    def __init__(self,
                 data,
                 yield_size:int,
                 
                 ):
        self.data = data
        self.yield_size = yield_size
        self._total_tasks = self.get_total_tasks()
    
    @abstractmethod
    async def generate_tasks(self) -> list[Task]:
        pass

    @abstractmethod
    def get_total_tasks(self) -> int:
        pass

class ProfileTaskGenerator(TaskGenerator):
    """This TaskGenerator generates FastProfileTask objects from a polars DataFrame. Each task profiles a BAM file."""
    def __init__(
        self,
        data: pl.LazyFrame,
        yield_size: int,
        container_engine: Engine,
        stb_file: str,
        profile_bed_file: str,
        gene_range_file: str,
        genome_length_file: str,
        num_procs: int = 4,
        polars_engine: str = "streaming",
    ) -> None:
        super().__init__(data, yield_size)
        self.stb_file = pathlib.Path(stb_file)
        self.profile_bed_file = pathlib.Path(profile_bed_file)
        self.gene_range_file = pathlib.Path(gene_range_file)
        self.genome_length_file = pathlib.Path(genome_length_file)
        self.num_procs = num_procs
        self.engine = container_engine
        self.polars_engine = polars_engine
        if type(self.data) is not pl.LazyFrame:
            raise ValueError("data must be a polars LazyFrame.")
        for path_attr in [
            self.stb_file,
            self.profile_bed_file,
            self.gene_range_file,
            self.genome_length_file,
        ]:
            if not path_attr.exists():
                raise FileNotFoundError(f"File {path_attr} does not exist.")

    def get_total_tasks(self) -> int:
        """Returns total number of profiles to be generated."""
        return self.data.select(size=pl.len()).collect(engine="streaming")["size"][0]
    
    async def generate_tasks(self) -> list[Task]:
        """Yeilds lists of FastProfileTask objects based on the data in batches of yield_size. This method yields the control back to the event loop
        while polars is collecting data to avoid blocking.
        """
        for offset in range(0, self._total_tasks, self.yield_size):
            batch_df = await self.data.slice(offset, self.yield_size).collect_async(engine="streaming")
            tasks = []
            for row in batch_df.iter_rows(named=True):
                inputs = {
                "bam-file": FileInput(row["bamfile"]),
                "sample-name": StringInput(row["sample_name"]),
                "stb-file": FileInput(self.stb_file),
                "bed-file": FileInput(self.profile_bed_file),
                "gene-range-table": FileInput(self.gene_range_file),
                "genome-length-file": FileInput(self.genome_length_file),
                "num-threads": IntInput(self.num_procs),
                }
                expected_outputs ={
                "profile":  FileOutput(row["sample_name"]+".parquet" ),
                "breadth":  FileOutput(row["sample_name"]+".breadth.parquet" ),
                }
                task = ProfileBamTask(id=row["sample_name"], inputs=inputs, expected_outputs=expected_outputs, engine=self.engine)
                tasks.append(task)
            yield tasks

class CompareTaskGenerator(TaskGenerator):
    """This TaskGenerator generates FastCompareTask objects from a polars DataFrame. Each task compares two profiles using compare_genomes functionality in
    zipstrain.compare module.
    
    Args:
        data (pl.LazyFrame): Polars LazyFrame containing the data for generating tasks.
        yield_size (int): Number of tasks to yield at a time.
        comp_config (database.GenomeComparisonConfig): Configuration for genome comparison.
        memory_mode (str): Memory mode for the comparison task. Default is "heavy".
        polars_engine (str): Polars engine to use. Default is "streaming".
        chrom_batch_size (int): Chromosome batch size for the comparison task in light memory mode. Default is 10000.
    """
    def __init__(
        self,
        data: pl.LazyFrame,
        yield_size: int,
        container_engine: Engine,
        comp_config: database.GenomeComparisonConfig,
        memory_mode: str = "heavy",
        polars_engine: str = "streaming",
        chrom_batch_size: int = 10000,
    ) -> None:
        super().__init__(data, yield_size)
        self.comp_config = comp_config
        self.engine = container_engine
        self.memory_mode = memory_mode
        self.polars_engine = polars_engine
        self.chrom_batch_size = chrom_batch_size
        if type(self.data) is not pl.LazyFrame:
            raise ValueError("data must be a polars LazyFrame.")
        
    def get_total_tasks(self) -> int:
        """Returns total number of pairwise comparisons to be made."""
        return self.data.select(size=pl.len()).collect(engine="streaming")["size"][0]

    async def generate_tasks(self) -> list[Task]:
        """Yeilds lists of FastCompareTask objects based on the data in batches of yield_size. This method yields the control back to the event loop
        while polars is collecting data to avoid blocking.
        """
        for offset in range(0, self._total_tasks, self.yield_size):
            batch_df = await self.data.slice(offset, self.yield_size).collect_async(engine="streaming")
            tasks = []
            for row in batch_df.iter_rows(named=True):
                inputs = {
                "mpile_1_file": FileInput(row["profile_location_1"]),
                "mpile_2_file": FileInput(row["profile_location_2"]),
                "scaffold_1_file": FileInput(row["scaffold_location_1"]),
                "scaffold_2_file": FileInput(row["scaffold_location_2"]),
                "null_model_file": FileInput(self.comp_config.null_model_loc),
                "stb_file": FileInput(self.comp_config.stb_file_loc),
                "min_cov": IntInput(self.comp_config.min_cov),
                "min-gene-compare-len": IntInput(self.comp_config.min_gene_compare_len),
                "memory-mode": StringInput(self.memory_mode),
                "chrom-batch-size": IntInput(self.chrom_batch_size),
                "genome-name": StringInput(self.comp_config.scope),
                "engine": StringInput(self.polars_engine),
                }
                expected_outputs ={
                "output-file":  FileOutput(row["sample_name_1"]+"_"+row["sample_name_2"]+"_comparison.parquet" ),

                }
                task = FastCompareTask(id=row["sample_name_1"]+"_"+row["sample_name_2"], inputs=inputs, expected_outputs=expected_outputs, engine=self.engine)
                tasks.append(task)
            yield tasks


    

class MapReadsTaskGenerator(TaskGenerator):
    pass

class Batch(ABC):
    """Batch is a collection of tasks to be executed as a group. This is a base class and should not be instantiated directly.
    A batch is the unit of execution meaning that the enitre batch is either run locally or submitted to a job scheduler like Slurm.
    
    Args:
        tasks (list[Task]): List of Task objects to be included in the batch.
        id (str): Unique identifier for the batch.
        run_dir (pathlib.Path): Directory where the batch will be executed.
        expected_outputs (list[Output]): List of expected outputs for the batch.
    """
    TEMPLATE_CMD = ""

    def __init__(self, tasks: list[Task],
                 id: str,
                 run_dir: pathlib.Path,
                 expected_outputs: list[Output],
                 file_semaphore: asyncio.Semaphore| None = None
                 ) -> None:
        self.id = id
        self.tasks = tasks
        self.run_dir = pathlib.Path(run_dir)
        self.batch_dir = self.run_dir / self.id
        self.retry_count = 0
        self.expected_outputs = expected_outputs
        self.file_semaphore = file_semaphore
        for output in self.expected_outputs:
            if isinstance(output, BatchFileOutput):
                output.register_batch(self)
        self._status = self._get_initial_status()
        for task in self.tasks:
            task._batch_obj = self
            task.file_semaphore = self.file_semaphore
            for output in task.expected_outputs.values():
                output.register_task(task)
            task._status= task._get_initial_status()
            task.map_io()
        
        self._runner_obj:Runner = None
            

    
    def _get_initial_status(self) -> str:
        """Returns the initial status of the batch based on the presence of the batch directory."""
        if not self.batch_dir.exists():
            return Status.NOT_STARTED.value
        with open(self.batch_dir / ".status", mode="r") as f:
            status_as_written = f.read().strip()
        if status_as_written in (Status.DONE.value, Status.SUCCESS.value):
            all_ready = True
            try:
                for output in self.expected_outputs:
                    if not output.ready():
                        all_ready = False
                        break
            except Exception:
                all_ready = False

            if all_ready:
                return Status.SUCCESS.value
            else:
                return Status.FAILED.value

    def cleanup(self) -> None:
        """The base class defines if any cleanup is needed after batch success. By default, it does nothing."""
        return None

    @abstractmethod
    async def cancel(self) -> None:
        """Cancels the batch. This method should be implemented by subclasses."""
        ...
        
    def outputs_ready(self) -> bool:
        """Check if all BATCH-LEVEL expected outputs are ready."""
        try:
            for output in self.expected_outputs:
                if not output.ready():
                    return False
            return True
        except Exception:
            return False
        
    async def _collect_task_status(self) -> list[str]:
        """Collects the status of all tasks asynchronously."""
        return await asyncio.gather(*[task.get_status() for task in self.tasks])

    @abstractmethod
    async def run(self) -> None:
        """Runs the batch. This method should be implemented by subclasses."""
        ...

    @abstractmethod
    def _parse_job_id(self, sbatch_output: str) -> str:
        """Parses the job ID from the sbatch output. This method should be implemented by subclasses."""
        ...

    @property
    def status(self) -> str:
        """Returns the current status of the batch."""
        return self._status

    @property
    def stats(self) -> dict[str, str]:
        """Returns a dictionary of task IDs and their statuses."""
        return {task.id: task.status for task in self.tasks}
    
    async def update_status(self) -> str:
        """Updates the status of the batch by collecting the status of all tasks."""
        await self._collect_task_status()
    def _set_file_semaphore(self, file_semaphore: asyncio.Semaphore) -> None:
        self.file_semaphore = file_semaphore
        for task in self.tasks:
            task.file_semaphore = file_semaphore

class LocalBatch(Batch):
    """Batch that runs tasks locally in a single shell script."""
    TEMPLATE_CMD = "#!/bin/bash\n"

    def __init__(self, tasks, id, run_dir, expected_outputs) -> None:
        super().__init__(tasks, id, run_dir, expected_outputs)
        self._script = self.TEMPLATE_CMD + "\nset -o pipefail\n"
        self._proc: asyncio.subprocess.Process | None = None 


    async def run(self) -> None:
        """This method runs all tasks in the batch locally by creating a shell script and executing it."""
        if self.status != Status.SUCCESS and self.status != Status.FAILED.value:
            self.batch_dir.mkdir(parents=True, exist_ok=True)
            self._status = Status.RUNNING.value


            await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)

            for task in self.tasks:
                if task.status == Status.NOT_STARTED.value:
                    task.task_dir.mkdir(parents=True, exist_ok=True)  # Create task directory
                    await write_file(task.task_dir / ".status", Status.NOT_STARTED.value, self.file_semaphore)

            script_path = self.batch_dir / f"{self.id}.sh" # Path to the shell script for the batch

            script = self._script
            for task in self.tasks:
                if task.status == Status.NOT_STARTED.value or task.status == Status.FAILED.value:
                    script += f"\n{task.pre_run}\n{task.command}\n{task.post_run}\n"
            
            await write_file(script_path, script, self.file_semaphore)

            await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)

            self._proc = await asyncio.create_subprocess_exec(
                "bash", f"{self.id}.sh",
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=self.batch_dir,
            )
            try:
                out_bytes, err_bytes = await self._proc.communicate()
            except asyncio.CancelledError:
                if self._proc and self._proc.returncode is None:
                    self._proc.terminate()
                raise

            await write_file(self.batch_dir / f"{self.id}.out", out_bytes.decode(), self.file_semaphore)
            await write_file(self.batch_dir / f"{self.id}.err", err_bytes.decode(), self.file_semaphore)

            if self._proc.returncode == 0 and self.outputs_ready():
                self.cleanup()
                self._status = Status.SUCCESS.value
                await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)
            else:
                self._status = Status.FAILED.value
                await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)
        
        elif self.status == Status.SUCCESS.value and self.outputs_ready():
            self._status = Status.SUCCESS.value
        else:
            self._status = Status.FAILED.value
    
    def _parse_job_id(self, sbatch_output):
        return super()._parse_job_id(sbatch_output)

    def cleanup(self) -> None:
        super().cleanup()

    async def cancel(self) -> None:
        """Cancels the local batch by terminating the subprocess if it's running."""
        if self._proc and self._proc.returncode is None:
            self._proc.terminate()
            try:
                await asyncio.wait_for(self._proc.wait(), timeout=5.0)
            except asyncio.TimeoutError:
                self._proc.kill()
                await self._proc.wait()
            self._status = Status.FAILED.value
            await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)


class SlurmBatch(Batch):
    """Batch that submits tasks to a Slurm job scheduler.
    
     Args:
        tasks (list[Task]): List of Task objects to be included in the batch.
        id (str): Unique identifier for the batch.
        run_dir (pathlib.Path): Directory where the batch will be executed.
        expected_outputs (list[Output]): List of expected outputs for the batch.
        slurm_config (SlurmConfig): Configuration for Slurm job submission. Refer to SlurmConfig class for details."""
    TEMPLATE_CMD = "#!/bin/bash\n"

    def __init__(self, tasks, id, run_dir, expected_outputs, slurm_config: SlurmConfig) -> None:
        super().__init__(tasks, id, run_dir, expected_outputs)
        self._check_slurm_works()
        self.slurm_config = slurm_config
        self._script = self.TEMPLATE_CMD + self.slurm_config.to_slurm_args() + "\nset -o pipefail\n"
        self._job_id = None

    def _check_slurm_works(self) -> None:
        """Checks if Slurm commands are available on the system."""
        try:
            subprocess.run(["sbatch", "--version"], capture_output=True, text=True, check=True)
            subprocess.run(["sacct", "--version"], capture_output=True, text=True, check=True)
        except:
            raise EnvironmentError("Slurm does not seem to be available or configured properly on this system.")

    async def cancel(self) -> None:
        """Cancel a running or submitted Slurm job."""
        if self._job_id:
            proc = await asyncio.create_subprocess_exec(
                "scancel", self._job_id,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            await proc.wait()

        self._status = Status.FAILED.value
        await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)
        
    async def run(self) -> None:
        """This method submits the batch to Slurm by creating a batch script and using sbatch command. It also monitors the job status until completion.
        This method is unavoidably different from LocalBatch.run() because of the nature of Slurm job submission.
        """
        
        if self.status != Status.SUCCESS and self.status != Status.FAILED.value:
            self.batch_dir.mkdir(parents=True, exist_ok=True)
            self._status = Status.RUNNING.value
            await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)
            # create task directories and initialize .status if needed
            
            for task in self.tasks:
                if task.status == Status.NOT_STARTED.value:
                    task.task_dir.mkdir(parents=True, exist_ok=True)
                    await write_file(task.task_dir / ".status", Status.NOT_STARTED.value, self.file_semaphore)
            # write the batch script (all tasks included)
            
            batch_path = self.batch_dir / f"{self.id}.batch"
            script=self._script
            for task in self.tasks:
                if task.status == Status.NOT_STARTED.value:
                    script += f"\n{task.pre_run}\n{task.command}\n{task.post_run}\n"
            
            await write_file(batch_path, script, self.file_semaphore)

            proc = await asyncio.create_subprocess_exec(
                "sbatch","--parsable", batch_path.name,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=str(self.batch_dir),
            )
            out_bytes, out_err = await proc.communicate()
            out = out_bytes.decode().strip() if out_bytes else ""
            if proc.returncode == 0:
                try:
                    self._job_id = self._parse_job_id(out)
                    self._status = Status.SUBMITTED.value
                    await self._wait_to_finish()
                except Exception:
                    self._status = Status.FAILED.value
            else:
                self._status = Status.FAILED.value
            
            if self._status == Status.SUCCESS.value and self.outputs_ready():
                self.cleanup()
                self._status = Status.SUCCESS.value
                await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)
            else:
                self._status = Status.FAILED.value
                await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)
            
        else:
            if self.status == Status.SUCCESS.value and self.outputs_ready():
                self._status = Status.SUCCESS.value
            else:
                self._status = Status.FAILED.value
            
    def _parse_job_id(self, sbatch_output: str) -> str:
        if match := re.search(r"(\d+)", sbatch_output):
            return match.group(1)
        else:
            raise ValueError("Could not parse job ID from sbatch output.")

    async def _wait_to_finish(self,sleep_duration:float=1.0):
        while self.status not in (Status.SUCCESS.value, Status.FAILED.value):
            await self.update_status()
            await asyncio.sleep(sleep_duration)
        
    async def update_status(self):
        if self._job_id is None:
            self._status=Status.NOT_STARTED.value
        else:
            await self._collect_task_status()
            out= await asyncio.create_subprocess_exec(
                "sacct", "-j", self._job_id, "--format=State", "--noheader","--allocations",
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            out_bytes, _ = await out.communicate()
            if out_bytes:
                state = out_bytes.decode().strip()
                if state in ["FAILED", "CANCELLED", "TIMEOUT"]:
                    self._status = Status.FAILED.value
                elif state=="RUNNING":
                    self._status = Status.RUNNING.value
                elif state in ["COMPLETED", "COMPLETING"]:
                    self._status = Status.SUCCESS.value
                else:
                    self._status = Status.PENDING.value
                await write_file(self.batch_dir / ".status", self._status, self.file_semaphore)

console = Console()

def get_cpu_usage():
    """Returns the current CPU usage percentage."""
    return psutil.cpu_percent(interval=0.1)

def get_memory_usage():
    """Returns the current memory usage percentage."""
    return psutil.virtual_memory().percent


class Runner(ABC):
    """Base Runner class to manage task generation, batching, and execution.
    
    Args:
        run_dir (str | pathlib.Path): Directory where the runner will operate.
        task_generator (TaskGenerator): An instance of TaskGenerator to produce tasks.
        container_engine (Engine): An instance of Engine to wrap task commands.
        batch_factory (Batch): The class that creates Batch instances. It should be a subclass of Batch or its subclasses.
        final_batch_factory (Batch): A callable that creates the final Batch instance.
        max_concurrent_batches (int): Maximum number of batches to run concurrently. Default is 1.
        poll_interval (float): Time interval in seconds to poll for batch status updates. Default is 1.0.
        tasks_per_batch (int): Number of tasks to include in each batch. Default is
        batch_type (str): Type of batch to use ("local" or "slurm"). Default is "local".
        slurm_config (SlurmConfig | None): Configuration for Slurm batches if batch_type
            is "slurm". Default is None.
        
    """
    TERMINAL_BATCH_STATES = {Status.SUCCESS.value, Status.FAILED.value}
    def __init__(self,
                    run_dir: str | pathlib.Path,
                    
                    task_generator: TaskGenerator,
                    container_engine: Engine,
                    batch_factory: Batch,
                    final_batch_factory: Batch,
                    max_concurrent_batches: int = 1,
                    poll_interval: float = 1.0,
                    tasks_per_batch: int = 10,
                    batch_type: str = "local",
                    slurm_config: SlurmConfig | None = None,
                    max_retries: int = 3,
                    ) -> None:
        self.run_dir = pathlib.Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.task_generator = task_generator
        self.container_engine = container_engine
        self.max_concurrent_batches = max_concurrent_batches
        self.poll_interval = poll_interval
        self.tasks_per_batch = tasks_per_batch
        self.batch_type = batch_type
        self.slurm_config = slurm_config
        self.tasks_queue: asyncio.Queue = asyncio.Queue(maxsize=2 * max_concurrent_batches * tasks_per_batch)
        self.batches_queue: asyncio.Queue = asyncio.Queue(maxsize=2 * max_concurrent_batches)
        self._finished_batches_count = 0
        self._success_batches_count = 0
        self._produced_tasks_count = 0
        self._active_batches: list[Batch] = []
        self._batch_counter = 0
        self._batcher_done = False
        self._final_batch_created = False
        self.batch_factory = batch_factory
        self.final_batch_factory = final_batch_factory
        self._failed_batches_count = 0
        self.max_retries = max_retries
        self.total_expected_tasks = self.task_generator.get_total_tasks()
        self.total_expected_batches = (self.total_expected_tasks + tasks_per_batch - 1) // tasks_per_batch
        self._shutdown_event = asyncio.Event()
        self._shutdown_initiated = False
    
    async def _refill_tasks(self):
        """Repeatedly call task_generator until it returns an empty list. This feeds tasks into the tasks_queue and waits for the queue to have space if it's full in an async manner."""
        async for tasks in self.task_generator.generate_tasks():
            for task in tasks:
                await self.tasks_queue.put(task)
                self._produced_tasks_count += 1
        await self.tasks_queue.put(None)
        
    @abstractmethod
    async def _batcher(self):
        ...

    async def _shutdown(self):
        """Cancel all active batches and signal shutdown."""
        if self._shutdown_initiated:
            return
        self._shutdown_initiated = True
        console.print("[yellow]Shutdown requested. Cancelling active jobs...[/]")

        for batch in list(self._active_batches):
            try:
                await batch.cancel()  
            except Exception as e:
                console.print(f"[red]Error cancelling batch {batch.id}: {e}[/]")

        # Signal the main loop to stop
        self._shutdown_event.set()

    async def run(self):
        """
        Run the producer, batcher and worker coroutines and present a live UI while working.
        Runs the task generator to produce tasks, batches them using the batcher,
        and executes batches with up to [max_concurrent_batches] parallel workers.
        UI: displays an overall panel (produced/finished counts), active batch Progress bars,
        and system stats (CPU/RAM) using Rich Live to mirror the Runner presentation.

        """
        asyncio.create_task(self._batcher())
        asyncio.create_task(self._refill_tasks())
        semaphore = asyncio.Semaphore(self.max_concurrent_batches)
        file_semaphore = asyncio.Semaphore(20)
        async def run_batch(batch: Batch):
            async with semaphore:
                while batch.status != Status.SUCCESS.value and batch.retry_count < self.max_retries:
                    await batch.run()
                    if batch.status == Status.SUCCESS.value:
                        break
                    else:
                        batch.retry_count += 1
                
                self._finished_batches_count += 1
                
                if batch.status == Status.SUCCESS.value:
                    self._success_batches_count += 1
                
                elif batch.status == Status.FAILED.value:
                    self._failed_batches_count += 1

                if batch in self._active_batches:
                    self._active_batches.remove(batch)

        # Rich progress objects
        
        overall_progress = Progress(
            TextColumn(f"[bold white]{type(self).__name__}[/]"),
            BarColumn(),
            TextColumn("• {task.fields[produced_tasks]}/{task.fields[total_expected_tasks]} tasks produced"),
            TextColumn("• {task.fields[finished_batches]}/{task.fields[total_expected_batches]} batches finished • {task.fields[failed_batches]} batches failed"),
            TimeElapsedColumn(),
            expand=True,
        )
        overall_task = overall_progress.add_task("overall", produced_tasks=0, total_expected_tasks=self.total_expected_tasks, finished_batches=0, total_expected_batches=self.total_expected_batches, failed_batches=0)

        batch_progress = Progress(
            TextColumn("[bold white]{task.fields[batch_id]}[/]"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("• {task.fields[status]}"),
            TimeElapsedColumn(),
            expand=True,
        )

        batch_to_progress_id: dict[Batch, int] = {}
        batch_task_totals: dict[Batch, int] = {}

        body = Panel(Group(
            Align.center(f"[bold magenta]ZipStrain {type(self).__name__}[/]\n", vertical="middle"),
            Panel(overall_progress, title="Overall Progress"),
            Panel(batch_progress, title="Active Batches", height=10),
            Panel(self._make_system_stats_panel(), title="System Stats", expand=True),
        ), border_style="magenta")
        
        loop = asyncio.get_running_loop()
        for sig in (signal.SIGINT, signal.SIGTERM):
            loop.add_signal_handler(sig, lambda s=sig: asyncio.create_task(self._shutdown()))
            
        with Live(body, console=console, refresh_per_second=2) as live:
            while not self._shutdown_event.is_set():
                await self._update_statuses()
                if self._batcher_done and self.batches_queue.empty() and len(self._active_batches) == 0:
                    if not self._final_batch_created:
                        final_batch = self._create_final_batch()
                        await self.batches_queue.put(final_batch)
                        self._final_batch_created = True
                    else:
                        break
                while len(self._active_batches) < self.max_concurrent_batches and not self.batches_queue.empty():
                    batch = await self.batches_queue.get()
                    if batch is not None:
                        batch._set_file_semaphore(file_semaphore)
                        self._active_batches.append(batch)
                        asyncio.create_task(run_batch(batch))
                # Update overall progress fields
                overall_progress.update(overall_task, produced_tasks=self._produced_tasks_count, finished_batches=self._finished_batches_count, failed_batches=self._failed_batches_count)  
                # Add newly queued batches into UI
                for batch in list(self._active_batches):
                    if batch not in batch_to_progress_id and batch.status not in self.TERMINAL_BATCH_STATES:
                        total = len(batch.tasks) if batch.tasks else 1
                        task_id = batch_progress.add_task("", total=total, batch_id=f"Batch {batch.id}", status=batch.status)
                        batch_to_progress_id[batch] = task_id
                        batch_task_totals[batch] = total
                # Remove finished batches from UI
                for batch, tid in list(batch_to_progress_id.items()):
                    if batch.status in self.TERMINAL_BATCH_STATES:
                        try:
                            batch_progress.remove_task(tid)
                        except Exception:
                            pass
                        del batch_to_progress_id[batch]
                        if batch in batch_task_totals:
                            del batch_task_totals[batch]
                # Update per-batch progress
                for batch, tid in batch_to_progress_id.items():
                    completed = sum(1 for t in batch.tasks if t.status in self.TERMINAL_BATCH_STATES)
                    total = batch_task_totals.get(batch, max(1, len(batch.tasks)))
                    batch_progress.update(tid, completed=completed, total=total, status=batch.status)
                # Update system panel
                body = Panel(Group(
                    Align.center(f"[bold magenta]ZipStrain {type(self).__name__}[/]\n", vertical="middle"),
                    Panel(overall_progress, title="Overall Progress"),
                    Panel(batch_progress, title="Active Batches"),
                    Panel(self._make_system_stats_panel(), title="System Stats", expand=True),
                ), border_style="magenta")
                live.update(body)
                await asyncio.sleep(self.poll_interval)

        # final UI summary
        console.clear()
        total_batches = self._batch_counter + (1 if self._final_batch_created else 0)
        summary = Panel(
            f"[bold green]Run finished![/]\n\n{self._success_batches_count}/{total_batches} batches succeeded.\n\nProduced tasks: {self._produced_tasks_count}\nElapsed: (see time in UI)",
            expand=True,
            title="Summary",
            border_style="green",
        )
        console.print(summary)
    
    async def _update_statuses(self):
        await asyncio.gather(*[batch.update_status() for batch in self._active_batches if batch.status not in self.TERMINAL_BATCH_STATES])

class CompareRunner(Runner):
    """
    Creates and schedules batches of FastCompareTask tasks using either local or Slurm batches.
    
    Args:
        run_dir (str | pathlib.Path): Directory where the runner will operate.
        task_generator (TaskGenerator): An instance of TaskGenerator to produce tasks.
        container_engine (Engine): An instance of Engine to wrap task commands.
        max_concurrent_batches (int): Maximum number of batches to run concurrently. Default is 1.
        poll_interval (float): Time interval in seconds to poll for batch status updates. Default is 1.0.
        tasks_per_batch (int): Number of tasks to include in each batch. Default is 10.
        batch_type (str): Type of batch to use ("local" or "slurm"). Default is "local".
        slurm_config (SlurmConfig | None): Configuration for Slurm batches if batch_type
            is "slurm". Default is None.    
    """

    def __init__(
        self,
        run_dir: str | pathlib.Path,
        task_generator: TaskGenerator,
        container_engine: Engine,
        max_concurrent_batches: int = 1,
        poll_interval: float = 1.0,
        tasks_per_batch: int = 10,
        batch_type: str = "local",
        slurm_config: SlurmConfig | None = None,
    ) -> None:
        if batch_type == "slurm":
            if slurm_config is None:
                raise ValueError("Slurm config must be provided for slurm batch type.")
            batch_factory = FastCompareSlurmBatch
            final_batch_factory = PrepareCompareGenomeRunOutputsSlurmBatch
        else:
            batch_factory = FastCompareLocalBatch
            final_batch_factory = PrepareCompareGenomeRunOutputsLocalBatch
        super().__init__(
            run_dir=run_dir,
            task_generator=task_generator,
            container_engine=container_engine,
            batch_factory=batch_factory,
            final_batch_factory=final_batch_factory,
            max_concurrent_batches=max_concurrent_batches,
            poll_interval=poll_interval,
            tasks_per_batch=tasks_per_batch,
            batch_type=batch_type,
            slurm_config=slurm_config,
        )




    async def _batcher(self):
        """
        Defines the batcher coroutine that collects tasks from the tasks_queue, groups them into batches,
        and puts the batches into the batches_queue. Each batch includes a CollectComps task to merge the outputs of the tasks in the batch.
        """
        buffer: list[Task] = []
        while True:
            task = await self.tasks_queue.get()
            if task is None:
                if buffer:
                    batch_id = f"batch_{self._batch_counter}"
                    self._batch_counter += 1
                    batch_tasks = buffer + [
                        CollectComps(
                            "concat_parquet",
                            {},
                            {"output-file": FileOutput(f"Merged_batch_{batch_id}.parquet")},
                            engine=self.container_engine,
                        )
                    ]
                    expected_outputs = [BatchFileOutput(f"concat_parquet/Merged_batch_{batch_id}.parquet")]
                    if self.batch_type == "slurm":
                        batch = self.batch_factory(
                            tasks=batch_tasks,
                            id=batch_id,
                            run_dir=self.run_dir,
                            expected_outputs=expected_outputs,
                            slurm_config=self.slurm_config,
                        )
                    else:
                        batch = self.batch_factory(
                            tasks=batch_tasks,
                            id=batch_id,
                            run_dir=self.run_dir,
                            expected_outputs=expected_outputs,
                        )
                    await self.batches_queue.put(batch)
                await self.batches_queue.put(None)
                self._batcher_done = True
                break
            buffer.append(task)
            if len(buffer) == self.tasks_per_batch:
                batch_id = f"batch_{self._batch_counter}"
                self._batch_counter += 1
                batch_tasks = buffer + [
                    CollectComps(
                        "concat_parquet",
                        {},
                        {"output-file": FileOutput(f"Merged_batch_{batch_id}.parquet")},
                        engine=self.container_engine,
                    )
                ]
                expected_outputs = [BatchFileOutput(f"concat_parquet/Merged_batch_{batch_id}.parquet")]
                if self.batch_type == "slurm":
                    batch = self.batch_factory(
                        tasks=batch_tasks,
                        id=batch_id,
                        run_dir=self.run_dir,
                        expected_outputs=expected_outputs,
                        slurm_config=self.slurm_config,
                    )
                else:
                    batch = self.batch_factory(
                        tasks=batch_tasks,
                        id=batch_id,
                        run_dir=self.run_dir,
                        expected_outputs=expected_outputs,
                    )
                await self.batches_queue.put(batch)
                buffer = []



    def _make_system_stats_panel(self):
        """helpers to create a system stats panel for the live UI."""
        def usage_bar(label: str, percent: float, color: str):
            p = Progress(
                TextColumn(f"[bold]{label}[/]"),
                BarColumn(bar_width=None, complete_style=color),
                TextColumn(f"{percent:.1f}%"),
                expand=True,
            )
            p.add_task("", total=100, completed=percent)
            return Panel(p, expand=True, width=30)

        cpu = psutil.cpu_percent(interval=None)
        ram = psutil.virtual_memory().percent
        cpu_panel = usage_bar("CPU", cpu, "cyan")
        ram_panel = usage_bar("RAM", ram, "magenta")
        return Columns([cpu_panel, ram_panel], expand=True, equal=True, align="center")

    def _create_final_batch(self) -> Batch:
        """Creates the final batch that prepares the overall outputs after all comparison batches are done."""
        final_task = PrepareCompareGenomeRunOutputs(
            id="prepare_outputs",
            inputs={"output-dir": StringInput("Outputs")},
            expected_outputs={},
            engine=self.container_engine,
        )
        expected_outputs = [BatchFileOutput("all_comparisons.parquet")]
        if self.batch_type == "slurm":
            final_batch=self.final_batch_factory(
                tasks=[final_task],
                id="Outputs",
                run_dir=self.run_dir,
                expected_outputs=expected_outputs,
                slurm_config=self.slurm_config,
            )
            final_batch._runner_obj = self
            return final_batch
        else:
            final_batch = self.final_batch_factory(
                tasks=[final_task],
                id="Outputs",
                run_dir=self.run_dir,
                expected_outputs=expected_outputs,
            )
            final_batch._runner_obj = self
            return final_batch
            


class ProfileBamTask(Task):
    """A Task that generates a mpileup file and genome breadth file in parquet format for a given BAM file using the fast_profile profile_bam command.
    The inputs to this task includes:

        - bam-file: The input BAM file to be profiled.

        - bed-file: The BED file specifying the regions to profile.
        
        - sample-name: The name of the sample being processed.

        - gene-range-table: A BED file specifying the gene ranges for the sample.

        - num-threads: The number of threads to use for processing.

        - genome-length-file: A file containing the lengths of the genomes in the reference fasta.

        - stb-file: The STB file used for profiling.

    Args:
        id (str): Unique identifier for the task.
        inputs (dict[str, Input]): Dictionary of input parameters for the task.
        expected_outputs (dict[str, Output]): Dictionary of expected outputs for the task.
        engine (Engine): Container engine to wrap the command."""
    
    TEMPLATE_CMD="""
    samtools index <bam-file>
    fastprofiler.sh <bed-file> <bam-file> <sample-name> <gene-range-table> <num-threads> 
    samtools idxstats <bam-file> |  awk '$3 > 0 {print $1}' > <sample-name>.parquet.scaffolds
    zipstrain utilities genome_breadth_matrix --profile <sample-name>.parquet \
        --genome-length <genome-length-file> \
        --stb <stb-file> \
        --min-cov <breadth-min-cov> \
        --output-file <sample-name>_breadth.parquet
    """
    
class FastCompareTask(Task):
    """A Task that performs a fast genome comparison using the fast_profile compare single_compare_genome command.
    
    Args:
        id (str): Unique identifier for the task.
        inputs (dict[str, Input]): Dictionary of input parameters for the task.
        expected_outputs (dict[str, Output]): Dictionary of expected outputs for the task.
        engine (Engine): Container engine to wrap the command.
        """
    TEMPLATE_CMD="""
    zipstrain compare single_compare_genome --mpileup-contig-1 <mpile_1_file> \
    --mpileup-contig-2 <mpile_2_file> \
    --scaffolds-1 <scaffold_1_file> \
    --scaffolds-2 <scaffold_2_file> \
    --null-model <null_model_file> \
    --stb-file <stb_file> \
    --min-cov <min_cov> \
    --min-gene-compare-len <min-gene-compare-len> \
    --memory-mode <memory-mode> \
    --chrom-batch-size <chrom-batch-size> \
    --output-file <output-file> \
    --genome <genome-name> \
    --engine <engine> 
    """


class CollectComps(Task):
    """A Task that collects and merges comparison parquet files from multiple FastCompareTask tasks into a single parquet file.
    
    Args:
        id (str): Unique identifier for the task.
        inputs (dict[str, Input]): Dictionary of input parameters for the task.
        expected_outputs (dict[str, Output]): Dictionary of expected outputs for the task.
        engine (Engine): Container engine to wrap the command."""
    TEMPLATE_CMD="""
    mkdir -p comps
    cp */*_comparison.parquet comps/
    zipstrain utilities merge_parquet --input-dir comps --output-file <output-file>
    rm -rf comps
    """
    
    @property
    def pre_run(self) -> str:
        return f"echo {Status.RUNNING.value} > {self.task_dir.absolute()}/.status"
    


class PrepareCompareGenomeRunOutputs(Task):
    """A Task that prepares the final output by merging all parquet files after all genome comparisons are done."""
    TEMPLATE_CMD="""
    mkdir -p <output-dir>/comps
    find "$(pwd)" -type f -name "Merged_batch_*.parquet" -print0 | xargs -0 -I {} ln -s {} <output-dir>/comps/
    zipstrain utilities merge_parquet --input-dir <output-dir>/comps --output-file <output-dir>/all_comparisons.parquet
    rm -rf <output-dir>/comps
    """
    
    @property
    def pre_run(self) -> str:
        """Sets the task status to RUNNING and changes directory to the runner's run directory since this task may need to access multiple batch outputs."""
        return f"echo {Status.RUNNING.value} > {self.task_dir.absolute()}/.status && cd {self._batch_obj._runner_obj.run_dir.absolute()}"
    



class FastCompareLocalBatch(LocalBatch):
    """A LocalBatch that runs FastCompareTask tasks locally."""
    def cleanup(self) -> None:
        tasks_to_remove = [task for task in self.tasks if isinstance(task, FastCompareTask)]
        for task in tasks_to_remove:
            self.tasks.remove(task)
            shutil.rmtree(task.task_dir)

class FastCompareSlurmBatch(SlurmBatch):
    """A SlurmBatch that runs FastCompareTask tasks on a Slurm cluster. Maybe removed in future"""
    def cleanup(self) -> None:
        tasks_to_remove = [task for task in self.tasks if isinstance(task, FastCompareTask)]
        for task in tasks_to_remove:
            self.tasks.remove(task)
            shutil.rmtree(task.task_dir)

class PrepareCompareGenomeRunOutputsLocalBatch(LocalBatch):
    pass


class PrepareCompareGenomeRunOutputsSlurmBatch(SlurmBatch):
    pass

    
def lazy_run_compares(
    run_dir: str | pathlib.Path,
    container_engine: Engine,
    comps_db: database.GenomeComparisonDatabase|None = None,
    tasks_per_batch: int = 10,
    max_concurrent_batches: int = 1,
    poll_interval: float = 5.0,
    execution_mode: str = "local",
    slurm_config: SlurmConfig | None = None,
    memory_mode: str = "heavy",
    chrom_batch_size: int = 10000,
    polars_engine: str = "streaming"
) -> None:
    """A helper function to quickly set up and run a CompareRunner with given parameters.
    
    Args:
        run_dir (str | pathlib.Path): Directory where the runner will operate.
        container_engine (Engine): An instance of Engine to wrap task commands.
        comps_db (GenomeComparisonDatabase | None): An instance of GenomeComparisonDatabase containing comparison data.
        tasks_per_batch (int): Number of tasks to include in each batch. Default is 10.
        max_concurrent_batches (int): Maximum number of batches to run concurrently. Default is 1.
        poll_interval (float): Time interval in seconds to poll for batch status updates. Default is 5.0.
        execution_mode (str): Execution mode, either "local" or "slurm". Default is "local".
    """
    task_generator = CompareTaskGenerator(
        data=comps_db.to_complete_input_table(),
        yield_size=tasks_per_batch,
        container_engine=container_engine,
        comp_config=comps_db.config,
        memory_mode=memory_mode,
        polars_engine=polars_engine,
        chrom_batch_size=chrom_batch_size,
    )
    if execution_mode=="local":
        batch_type="local"
    elif execution_mode=="slurm":
        batch_type="slurm"
    else:
        raise ValueError(f"Unknown execution mode: {execution_mode}")
    runner = CompareRunner(
        run_dir=pathlib.Path(run_dir),
        task_generator=task_generator,
        container_engine=container_engine,
        max_concurrent_batches=max_concurrent_batches,
        poll_interval=poll_interval,
        tasks_per_batch=tasks_per_batch,
        batch_type=batch_type,
        slurm_config=slurm_config,
    )
    asyncio.run(runner.run())