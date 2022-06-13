from setuptools import setup, Extension
from Cython.Build import cythonize
import os
import numpy as np

if "HDF5_PATH" in os.environ:
    hdf5_dir = os.environ["HDF5_PATH"]
elif "YT_DEST" in os.environ:
    hdf5_dir = os.environ["YT_DEST"]
else:
    print("You need to set HDF5_PATH or YT_DEST in your environment.")
    raise RuntimeError
inc = os.path.join(hdf5_dir, "include")
lib = os.path.join(hdf5_dir, "lib")
ext = Extension(name = 'predator_prey_solver',
                 sources=[ "predator_prey_solver_run.pyx" ,
                    "predator_prey_solver.C", "BE_chem_solve.C"],
                 include_dirs=[inc, ".", np.get_include()],
                 library_dirs=[lib],
                 libraries=['m', 'hdf5', 'hdf5_hl', "stdc++"],
                 language="C", extra_compile_args=["-w", "-g"])
setup(
    ext_modules= cythonize([ext], language='C'),
)
