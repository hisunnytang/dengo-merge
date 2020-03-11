from evolve_free_fall import run_dengo_freefall, compare_dengo_grackle
from utilities import set_env_variables
import pytest
import os

#print("---------------------------")
#print(os.environ)
#print(os.environ["TRAVIS_BUILD_DIR"])
#print("--------------------------")

if "TRAVIS_BUILD_DIR" not in os.environ:
    set_env_variables("HDF5_DIR", "/home/kwoksun2/anaconda3")
    set_env_variables("CVODE_PATH", "/home/kwoksun2/cvode-3.1.0/instdir")
    set_env_variables("HDF5_PATH", "/home/kwoksun2/anaconda3")
    set_env_variables("SUITESPARSE_PATH", "/home/kwoksun2/SuiteSparse")
    set_env_variables("DENGO_INSTALL_PATH", "/home/kwoksun2/dengo_install")
else:
    # then we assume that the libraries are installed relative to the dengo
    # path, so the below paths are the relative default install path
    set_env_variables("HDF5_DIR", "hdf5_install")
    set_env_variables("CVODE_PATH", "cvode-3.1.0/instdir")
    set_env_variables("HDF5_PATH", "hdf5_install")
    set_env_variables("SUITESPARSE_PATH", "suitesparse")
    set_env_variables("DENGO_INSTALL_PATH", "dengo_install")


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

    if "TRAVIS_BUILD_DIR" not in os.environ:
        os.environ["DENGO_PATH"] = "/home/kwoksun2/dengo-git"

    # if you install it through the install scripts
    set_env_variables("HDF5_DIR", "hdf5_install")
    set_env_variables("CVODE_PATH", "cvode-3.1.0/instdir")
    set_env_variables("HDF5_PATH", "hdf5_install")
    set_env_variables("SUITESPARSE_PATH", "suitesparse")
    set_env_variables("DENGO_INSTALL_PATH", "dengo_install")


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
