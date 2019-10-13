from evolve_free_fall import run_dengo_freefall, compare_dengo_grackle
import pytest
import os


def set_env_variables(var, path):
    if var not in os.environ:
        os.environ[var] = path
        print("updating {} = {}".format(var, path))


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


def test_freefall(setup_solver_options):
    set_env_variables("HDF5_DIR", "/home/kwoksun2/anaconda3")
    set_env_variables("CVODE_PATH", "/home/kwoksun2/cvode-3.1.0/instdir")
    set_env_variables("HDF5_PATH", "/home/kwoksun2/anaconda3")
    set_env_variables("SUITESPARSE_PATH", "/home/kwoksun2/SuiteSparse")
    set_env_variables("DENGO_INSTALL_PATH", "/home/kwoksun2/dengo_install")
    run_dengo_freefall(setup_solver_options)
    compare_dengo_grackle(setup_solver_options["output_dir"])

def main():
   run_dengo_freefall(so)

if __name__ == "__main__":
    main()
