import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import sunny_solver_run
sunny_solver_run.main_run_sunny()