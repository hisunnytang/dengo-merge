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

`Dengo` can also be installed through `pypi` using `pip install dengo-chemistry`.
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



## Example
While we recommend our readers to read the [cookbook](https://hisunnytang.github.io/dengo-cookbook/landing.html) for a more detailed usage, here we provide a small and quick example with 1 reactions.
$$
\begin{split}
\begin{align*}
    \rm H + e^- &\rightarrow \rm H^+ + 2 e^- \quad &(k01)
\end{align*}
\end{split}
$$

### Defining Species
One would first need to define the chemical species using `AtomicSpecies`, and `MolecularSpecies`. The weight of the atom would be looked up from the periodic table internally in `dengo`. User has to specified the free electrons in each chemical species (hydrogen atom). `de` represents the free electrons.


```python
from dengo.reaction_classes import AtomicSpecies, MolecularSpecies
from dengo.chemical_network import species_registry

HI = AtomicSpecies('H', free_electrons=0)
HII = AtomicSpecies("H", free_electrons=1)
de  = species_registry['de']
```

The species would be automatically registered in `dengo`. and are accesible through `species_registry`.


### Defining Reactions

A reaction to dengo is primarily composed of:
- name of the reaction, i.e. k01, k02
- LHS of the reaction /input species
- RHS of the reaction / output species
- reaction rate, that often times dependent on the temperature of the gas parcel.

#### Reactants and Products

For example, dengo expects the LHS and RHS to have an input of the type for k01:
```
LHS = [(1, HI), (1, de)]
RHS = [(1, HII), (2, de)]
```

#### Reaction Rate
Reaction are oftentimes dependent on temperature. dengo expects a reaction rate function that take state as input. state contains not only temperature in $\rm K$, but also in Kelvin log scale, and in $\rm eV / K$  (electron volts per kelvin), and $\rm eV / K$ in log scale.

```python
from dengo.reaction_classes import reaction
from dengo.chemical_network import reaction_registry

tiny = 1e-60
# -- k01 --
@reaction('k01', [   (1,HI),   (1,de) ], [  (1,HII),   (2,de) ])
def rxn(state):
    vals = np.exp(-32.71396786375
                  + 13.53655609057*state.logtev
                  - 5.739328757388*state.logtev**2
                  + 1.563154982022*state.logtev**3
                  - 0.2877056004391*state.logtev**4
                  + 0.03482559773736999*state.logtev**5
                  - 0.00263197617559*state.logtev**6
                  + 0.0001119543953861*state.logtev**7
                  - 2.039149852002e-6*state.logtev**8)
    # taken from Abel 1999
    vals = np.maximum(vals , tiny *np.ones((len(state.T))) )
    return vals

```python
class state:
    def __init__(self, T_bounds=(1, 1e8), n_bins=1024):
        """Initialize the range of temperature over which the rate tables are generated

        Parameters
        ----------
        T_bounds: List[Float, Float], optional (default=(1,1e8))
            the range over which the rates table is interpolated
        n_bins: int, optional (default=1024)
        """""

        self.n_bins = n_bins
        self.T = np.logspace(
            np.log(T_bounds[0]), np.log(T_bounds[1]), n_bins, base=np.e

        )
        self.logT = np.log(self.T)
        self.tev = self.T / tevk
        self.logtev = np.log(self.tev)
        self.T_bounds = T_bounds
```

#### Creating `ChemicalNetwork`
We can now assemble them together in this ChemicalNetwork. This object helps us do all the neccessary computations to arrive at the symbolic rhs and jacobian functions, which ultimately eases us of the process of deriving them by hand.
```python
import dengo
from dengo.chemical_network import ChemicalNetwork

simpleNetwork = ChemicalNetwork()
simpleNetwork.add_reaction("k01")
simpleNetwork.init_temperature((1e0, 1e8))
```
This allows us to write the ODE solver from dengo templates

```python
solver_name = 'simpleNetwork'
output_dir  = '.'

simpleNetwork.write_solver(
    solver_name,
    output_dir=output_dir,
    solver_template="be_chem_solve/rates_and_rate_tables",
    ode_solver_source="BE_chem_solve.C",
    init_values=init_values,
)
```

```
simpleNetwork_solver.C
simpleNetwork_solver.h
simpleNetwork_solver_main.C
simpleNetwork_solver_main.py
simpleNetwork_solver_run.pxd
simpleNetwork_solver_run.pyx
simpleNetwork_solver_run.pyxbld
simpleNetwork_solver_run.pyxdep
simpleNetwork_tables.h5
```




### Compile the Cython Module

The solver templates come with a `.pyx` files that lets you build a python interface to the C-library. A more well-rounded Cython tutorial can be found here [Cython Tutorial](https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html)

```python
import pyximport
import numpy as np
solver_name = "simpleNetwork"

pyximport.install(
    setup_args={"include_dirs":np.get_include()},
    reload_support=True,
    inplace=True,
    language_level=3
)

simple_solver_run = pyximport.load_module(
    f"{solver_name}_solver_run",
    f"{solver_name}_solver_run.pyx",
    build_inplace = True,
    pyxbuild_dir = "_dengo_temp",
)
```

### Invoking the Cython solver
`{solver_name}_solver_run.run_{solver_name}(init_values, dt, niter)` is the entry point for the built solver. It takes an dictionary that contains the abundances, and thermal energy, and dt, time to advance the fluid parcel, niter the maximum number of iterations as arguments.

It assumes that abundances are in number density with the units of $\rm cm^{-3}$, and thermal energy density in $\rm erg/g$ (thermal energy per mass density) $\frac{1}{\rho} \sum_i \frac{n_i k T}{\gamma_i - 1}$. niter implicitly sets the initial timestep for the solver, i.e. $dt_{\rm solver} = \rm dt / \rm niter$ .

```python
NCELLS = 1
density = 1e-2

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values['H_1']     = init_array
init_values['H_2']     = init_array
init_values['de']      = init_array
init_values['ge']      = np.ones(NCELLS)*1e13

total_density = simpleNetwork.calculate_total_density(init_values)
init_values = simpleNetwork.convert_to_mass_density(init_values)
init_values['de'] = simpleNetwork.calculate_free_electrons(init_values)
init_values['density'] = simpleNetwork.calculate_total_density(init_values)
number_density = simpleNetwork.calculate_number_density(init_values)

rv, rv_int = simple_solver_run.run_simpleNetwork(init_values, 1e16, niter = 1e5)
```
