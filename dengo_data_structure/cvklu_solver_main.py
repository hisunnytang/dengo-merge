import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import cvklu_solver_run
cvklu_solver_run.main_run_cvklu()