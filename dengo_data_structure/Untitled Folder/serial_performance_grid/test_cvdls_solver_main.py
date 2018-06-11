import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import test_cvdls_solver_run
test_cvdls_solver_run.main_run_test_cvdls()