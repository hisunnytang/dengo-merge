import numpy
import pyximport
pyximport.install(setup_args={"include_dirs":numpy.get_include()},
                  reload_support=True)

import cvspils_9species_cooling_solver_run
cvspils_9species_cooling_solver_run.main_run_cvspils_9species_cooling()