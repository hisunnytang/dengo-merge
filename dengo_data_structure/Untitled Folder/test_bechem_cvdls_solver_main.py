import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import test_bechem_cvdls_solver_run
test_bechem_cvdls_solver_run.main_run_test_bechem_cvdls()