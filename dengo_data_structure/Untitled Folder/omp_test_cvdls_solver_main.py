import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import omp_test_cvdls_solver_run
omp_test_cvdls_solver_run.main_run_omp_test_cvdls()