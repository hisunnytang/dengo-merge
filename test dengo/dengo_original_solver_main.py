import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import dengo_original_solver_run
dengo_original_solver_run.main_run_dengo_original()