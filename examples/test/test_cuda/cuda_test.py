from dengo.chemical_network import \
        ChemicalNetwork, \
        reaction_registry, \
        cooling_registry
import dengo.primordial_rates
import dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh, G
import os
import jinja2

os.environ['HDF5_PATH'] = '/home/kwoksun2/anaconda3'
os.environ['DENGO_INSTALL_PATH'] = '/home/kwoksun2/dengo_install'
os.environ['CVODE_PATH'] = "/home/kwoksun2/dengo-merge/cvode-3.1.0/instdir"
os.environ['SUITESPARSE_PATH'] = "/home/kwoksun2/SuiteSparse"


def setup_primordial_network():
    """Initial a ChemicalNetwork object
       for primordial network 9-species model
    Return:
        primordial: ChemicalNetwork with primordial reactions and cooling
    """
    # this register all the rates specified in `primordial_rates.py`
    dengo.primordial_rates.setup_primordial()

    # initialize the chmical network object
    primordial = ChemicalNetwork()

    # add all the reactions
    primordial.add_reaction("k01")
    primordial.add_reaction("k02")
    primordial.add_reaction("k03")
    primordial.add_reaction("k04")
    primordial.add_reaction("k05")
    primordial.add_reaction("k06")
    primordial.add_reaction("k07")
    primordial.add_reaction("k08")
    primordial.add_reaction("k09")
    primordial.add_reaction("k10")
    primordial.add_reaction("k11")
    primordial.add_reaction("k12")
    primordial.add_reaction("k13")
    primordial.add_reaction("k14")
    primordial.add_reaction("k15")
    primordial.add_reaction("k16")
    primordial.add_reaction("k17")
    primordial.add_reaction("k18")
    primordial.add_reaction("k19")
    primordial.add_reaction("k21")
    primordial.add_reaction("k22")
    primordial.add_reaction("k23")

    primordial.add_cooling("brem")
    primordial.add_cooling("reHII")
    primordial.add_cooling("reHeIII")
    primordial.add_cooling("gloverabel08")
    primordial.add_cooling("ceHI")
    primordial.add_cooling("h2formation")
    primordial.add_cooling("reHeII2")
    primordial.add_cooling("reHeII1")
    primordial.add_cooling("ciHeIS")
    primordial.add_cooling("ceHeII")
    primordial.add_cooling("ciHI")
    primordial.add_cooling("ceHeI")
    primordial.add_cooling("gammah")
    primordial.add_cooling("ciHeI")
    primordial.add_cooling("ciHeII")
    primordial.add_cooling("cie_cooling")
    primordial.add_cooling("compton")

    # This defines the temperature range for the rate tables
    primordial.init_temperature((1e0, 1e8))

    #primordial.enforce_conservation = True
    #primordial.set_equilibrium_species("H2_2")

    return primordial


def write_network(network):

    network.write_solver("primordial",
            solver_template="cv_omp/sundials_CVDls",
            ode_solver_source="initialize_cvode_solver.C",
            output_dir = 'suitesparse')

def run_suitesparse_test():
    os.chdir("suitesparse")
    os.system("make -DNTHREADS=24")

    h5path = os.environ['HDF5_PATH']
    dengopath = os.environ['DENGO_INSTALL_PATH']
    cvodepath = os.environ['CVODE_PATH']

    os.system(f"{CC} test_performance.C -I{dengopath}/include -I{cvodepath}/include -L{dengopath} -ldengo -L{h5path} -lhdf5 -lhdf5_hl -lstdc++ -lm")

if __name__ == '__main__':
    network = setup_primordial_network()
    network.write_cuda_solver(solver_name="primordial_cuda", output_dir='cuda_solver')
    write_network(network)
