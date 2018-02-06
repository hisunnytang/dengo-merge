import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import bechem_9species_solver_run
bechem_9species_solver_run.main_run_bechem_9species()