from evolve_free_fall import run_dengo_freefall, compare_dengo_grackle
import pytest
import os


@pytest.fixture
def setup_solver_options(update_options={}):
    solver_options = {"output_dir": "temp_freefall",
                      "solver_name": "primordial",
                      "use_omp": True,
                      "use_cvode": True,
                      "use_suitesparse": True,
                      "niters": 1,
                      "NCELLS": 1,
                      "reltol": 1.0e-6}
    solver_options.update(update_options)
    return solver_options

@pytest.mark.skip(reason='seg fault on github action?')
def test_freefall(setup_solver_options):

    run_dengo_freefall(setup_solver_options)
    compare_dengo_grackle(setup_solver_options["output_dir"])

def main():
    solver_options = {"output_dir": "temp_freefall",
                      "solver_name": "primordial",
                      "use_omp": False,
                      "use_cvode": False,
                      "use_suitesparse": False,
                      "niters": 1,
                      "NCELLS": 1,
                      "reltol": 1.0e-6}
    run_dengo_freefall(solver_options)

if __name__ == "__main__":
    main()
