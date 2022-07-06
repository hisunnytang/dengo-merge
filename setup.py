import os
import os.path
import subprocess
import sys
import time

# read the contents of your README file
from pathlib import Path

import setuptools
from setuptools import find_packages, setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


VERSION = "0.1.1"


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )
    config.add_subpackage("dengo", "dengo")

    return config


setup(
    name="dengo",
    version=VERSION,
    description="A framework for creating chemical and radiative cooling"
    + "reaction networks",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: C",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="astronomy astrophysics chemistry",
    author="Sunny Tang, Matthew J. Turk, Devin Silvia",
    author_email="hisunnytang@gmail.com, matthewturk@gmail.com, devin.silvia@gmail.com",
    url="https://github.com/hisunnytang/dengo-merge",
    license="GPL-3",
    # configuration=configuration,
    zip_safe=False,
    install_requires=[
        "ChiantiPy>=0.12.0",
        "Cython>=0.29.30",
        "h5py>=3.6.0",
        "Jinja2>=3.0.3",
        "numpy>=1.18.5",
        "setuptools>=62.3.2",
        "sympy>=1.10.1",
        "matplotlib>=3.5.2",
        "pytest>=7.1.2",
        "unyt>=2.8.0",
        "docutils",
    ],
    packages=find_packages(include=["dengo", "dengo.*"]),
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    package_data={"": ["dengo/*.csv", "dengo/templates/*", "dengo/solvers/*"]},
)
