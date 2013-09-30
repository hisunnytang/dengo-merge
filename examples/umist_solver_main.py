import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import umist_solver_run
umist_solver_run.main_run_umist()