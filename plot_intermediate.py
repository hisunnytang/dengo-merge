import matplotlib.pyplot as plt
import numpy as np
import h5py
import glob
import argparse

class ResultsPlotter:

    def __init__(self, network_name):
        self.network_name = network_name
        self.grab_results()

    def grab_results(self):
        fns = glob.glob("%s_intermediate_*.h5" % (self.network_name))
        if len(fns) == 0: raise RuntimeError
        f = h5py.File(fns[0])
        self.species = f.keys()
        self.nspecies = len(self.species)
        self.ncells = f[self.species[0]].shape[0]
        self.data = dict(
            (s, np.empty((self.ncells, len(fns)), dtype="float64"))
            for s in self.species)
        f.close()
        self.t = np.empty(len(fns), dtype="float64")
        self.dt = np.empty(len(fns), dtype="float64")
        for i,fn in enumerate(sorted(fns)):
            if i % 1000 == 0:
                print "Loading from %s" % fn
            f = h5py.File(fn)
            for s in self.species:
                self.data[s][:, i] = f[s][:]
            self.t[i] = f["/"].attrs["t"]
            self.dt[i] = f["/"].attrs["timestep"]
            f.close()
        print "Loaded %s files" % (len(fns))

    def plot_all(self):
        plt.clf()
        for s in sorted(self.species):
            plt.loglog(self.t, self.data[species], '-', label=s)
        plt.legend(loc = "best")
        plt.savefig("%s_time.png" % (self.network_name))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("network_name", nargs=1)
    args = parser.parse_args()
    rp = ResultsPlotter(args.network_name[0])
    rp.plot_all()
