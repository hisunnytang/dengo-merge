import numpy
import pyximport

pyximport.install(setup_args={"include_dirs": numpy.get_include()}, reload_support=True)

import cvdls_9species_solver_run

cvdls_9species_solver_run.main_run_cvdls_9species()
