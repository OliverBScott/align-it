"""setup script for building + installing python wrappers."""
import pathlib
import sys
import os

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, pip 10 or greater is required,\n"
        " or you will need to install the PEP 518 requirements"
        " in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages
from distutils import sysconfig

# Find rdkit headers (assumes using a conda environment)
includes = pathlib.Path(sysconfig.get_python_inc()).parent
rdkit_include = includes / 'rdkit'

# Raise error if headers not found
if not rdkit_include.exists():
    raise RuntimeError(f"Could not find rdkit include headers at {rdkit_include},"
                       f" check that rdkit is installed correctly.")

# Get conda libs
conda_prefix = os.getenv("CONDA_PREFIX")
if not conda_prefix:
    conda_prefix = os.getenv("MINICONDAPATH")
if not conda_prefix:
    raise RuntimeError("Could not get environment variable CONDA_PREFIX")
conda_prefix = pathlib.Path(conda_prefix)
conda_libs = conda_prefix / 'lib'

# CMake's arguments (use rdkit and of course build wrappers)
CMAKE_ARGUMENTS = [
    "-DBUILD_RDKIT_SUPPORT=ON",
    "-DBUILD_PYTHON_SUPPORT=ON",
    f"-DRDKIT_INCLUDE_DIR={rdkit_include}/",
    f"-DCONDA_LIB_DIR={conda_libs}/",
]

setup(
    name="pyalignit",
    version="1.0.4",
    description="Python wrappers for the Align-it™ tool from Silicos-it",
    url="https://github.com/OliverBScott/align-it",
    author="OliverBScott",
    license="GPLv3",
    packages=find_packages(),
    cmake_install_dir="pyalignit",
    include_package_data=True,
    extras_require={"test": ["pytest"]},
    python_requires=">=3.6",
    cmake_args=CMAKE_ARGUMENTS,
    cmake_minimum_required_version="3.14"
)
