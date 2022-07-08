# Installation Guide
This chapter walks you through the steps needed to install `dengo`. `Dengo` can also be installed through `pypi` using `pip install dengo-chemistry`.

Although it is not neccessary, it is however highly recommended users to install it a seperate environment, such that it does not mess with your current python packages and installations.

While it can work standlone with the solver available in the shipped package, it is advisible to install `dengo` with the environement maintained by `conda`. This allows user to use the [Sundials CVODE](https://computing.llnl.gov/projects/sundials) with [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) sparse matrix algorithms. For dengo to function properly, environemtal variables have to be set accordingly for dengo to locate and link the require
packages `C` or `cython` modules.

```bash
# Clone the Current Repo
git clone https://github.com/hisunnytang/dengo-merge.git

# Create the Conda Environment
conda create -n dengo
# Activate the Environment
conda activate dengo

# Install all the neccessary packages
conda install -c conda-forge "sundials>=6.2"
conda install -c conda-forge ChiantiPy Cython h5py Jinja2 numpy setuptools sympy matplotlib pytest unyt pip

# Install dengo
cd dengo-merge
pip install -e .

# Define the environmental variables
export HDF5_PATH=${CONDA}/envs/test
export SUITESPARSE_PATH=${CONDA}/envs/test
export CVODE_PATH=${CONDA}/envs/test
export DENGO_INSTALL_PATH=${your_directory}/dengo_install
export LD_LIBRARY_PATH=${CONDA}/envs/test/lib:${DENGO_INSTALL_PATH}/lib:$LD_LIBRARY_PATH
```


[ChiantiPy](https://chiantipy.readthedocs.io/en/latest/) is a required package that contains the atomic database for astrophysical specroscopy data. Please read the their documentation or install-chiantipy.sh available in this repository.

[Sundials CVODE](https://computing.llnl.gov/projects/sundials/cvode) is a `LLNL` maintained solver for stiff and nonstiff ordinary differential equation (ODE) systems (initial value problem) given in explicit form $\frac{dy}{dt} = f(t,y)$. Current `dengo` templates work with CVODE with version greater thatn v6.0.0. If your reaction network has a sparse jacobian in nature, CVODE could take advantage of the sparsity, and use a sparse solver of Newton Iteration. [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) is a suite of sparse matrix algorithms that CVODE can take use to speed up the computation.

Note that CVODE and SuiteSparse are optional. While `Dengo` can function with the current installation, `Dengo` could make use of external ODE solver libraries to obtain additional performance gain.
