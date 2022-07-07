# Dengo: An Engineer for Chemistry Solvers

[![Latest Version](https://img.shields.io/pypi/v/dengo-chemistry?logo=dengo-chemistry)](https://pypi.org/project/dengo-chemistry/)
![Dengo Testing](https://github.com/hisunnytang/dengo-merge/actions/workflows/python-package.yml/badge.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)]

Hi there! Welcome to Dengo. Dengo is a Python system for symbolically describing a system of chemical species, chemical kinetic rate equations, cooling functions, and then producing from these a set of numerical kernels to calculate the “right hand side” and Jacobian of this system. These two numerical kernels can then be linked into one of several ODE solvers.

Dengo is best thought of as a way to create a chemistry and cooling solver that you can stick into a codebase. Rather than trying to tie together the four separate parts – species, reactions, rate coefficients and the solver itself – into a single whole, Dengo allows you to construct each item individually and join them at the final step before inserting them into a simulation code.

A online version of the cookbook and documentation can be found here
https://hisunnytang.github.io/dengo-cookbook/landing.html

For more information, please contact the authors:
- Sunny Tang   (hisunnytang@gmail.com)
- Matthew Turk (matthew.turk@gmail.com)
- Devin Silvia (devin.silvia@gmail.com)


## Installation
While dengo can work standlone with the solver available in the shipped package, it is advisible to install `dengo` with the environement maintained by `conda`. This allows user to use the [Sundials CVODE](https://computing.llnl.gov/projects/sundials) with [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) sparse matrix algorithms.

```bash
# Clone the Current Repo
git clone https://github.com/hisunnytang/dengo-merge.git

# Create the Conda Environment
conda create -n dengo
# Activate the Environment
conda activate dengo

# Install all the neccessary packages
conda install -c conda-forge "sundials>=6.2" ChiantiPy Cython h5py Jinja2 numpy setuptools sympy matplotlib pytest unyt pip

# Install dengo
cd dengo-merge
pip install -e .
```
