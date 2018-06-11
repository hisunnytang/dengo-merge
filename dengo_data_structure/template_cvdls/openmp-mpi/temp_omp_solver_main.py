import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import temp_omp_solver_run
temp_omp_solver_run.main_run_temp_omp()