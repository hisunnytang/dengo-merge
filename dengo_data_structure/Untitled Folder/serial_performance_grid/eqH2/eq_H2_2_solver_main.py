import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import eq_H2_2_solver_run
eq_H2_2_solver_run.main_run_eq_H2_2()