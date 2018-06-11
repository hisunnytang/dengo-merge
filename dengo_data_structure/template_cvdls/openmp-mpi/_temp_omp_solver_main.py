import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import _temp_omp_solver_run
_temp_omp_solver_run.main_run__temp_omp()