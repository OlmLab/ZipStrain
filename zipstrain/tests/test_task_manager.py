from zipstrain import task_manager
import pytest

def test_create_new_input_type():
    class NewInput(task_manager.Input):
        def validate(self) -> None:
            if not isinstance(self.value, str):
                raise ValueError(f"Input value {self.value} is not a string.")

        def get_value(self) -> str:
            return self.value

    with pytest.raises(ValueError, match="Input value 123 is not a string."):
        NewInput(123)  # This should raise a ValueError

def test_file_input(tmp_path):
    # Creation of an Input file that does not exist should raise FileNotFoundError
    with pytest.raises(FileNotFoundError, match="Input file .* does not exist."):
        task_manager.FileInput("non_existent_file.txt")
    
    #Now create a real file and test
    test_file = tmp_path / "test.txt"
    test_file.write_text("Just a test file.")
    file_input = task_manager.FileInput(str(test_file))
    assert file_input.get_value() == str(test_file.absolute())
    
def test_string_input():
    string_input = task_manager.StringInput("test_string")
    assert string_input.get_value() == "test_string"
    
    with pytest.raises(ValueError, match="Input value 123 is not a string."):
        task_manager.StringInput(123)  # This should raise a ValueError
        
def test_int_input():
    int_input = task_manager.IntInput(42)
    assert int_input.get_value() == "42"

    with pytest.raises(ValueError, match="Input value 'test_string' is not an integer."):
        task_manager.IntInput("test_string")  # This should raise a ValueError

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
     
    

