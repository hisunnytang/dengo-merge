import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import cache_solver_run
cache_solver_run.main_run_cache()