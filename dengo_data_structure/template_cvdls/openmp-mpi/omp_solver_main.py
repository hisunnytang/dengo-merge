import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import omp_solver_run
omp_solver_run.main_run_omp()