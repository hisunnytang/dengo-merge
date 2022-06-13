import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import predator_prey_solver_run
predator_prey_solver_run.main_run_predator_prey()