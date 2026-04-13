"""
ZipStrain
========================
ZipStrain is a bioinformatics tool for strain-level analysis of metagenomics data. It is designed with scalability and efficiency in mind, leveraging advanced data processing techniques to handle large datasets.

"""
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
import tomllib


def _resolve_version() -> str:
    try:
        return package_version("zipstrain")
    except PackageNotFoundError:
        pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
        with pyproject_path.open("rb") as handle:
            return tomllib.load(handle)["project"]["version"]


__version__ = _resolve_version()
