import numpy
import pyximport
import os

# write the solver for various network

# test the cython intallation
os.environ["HDF5_DIR"] = /home/kwoksun2/anaconda3
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import primordial_solver_run

# this runs the `primordial_main` function in the primordial_solver
# reading the primordial_initial_conditions.h5
primordial_solver_run.main_run_primordial()

# test the run_primordial(init, )