import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import h5py
import glob
import argparse
import docutils.utils.roman as roman

class ResultsPlotter:

    def __init__(self, network_name):
        self.network_name = network_name
        self.grab_results()

    def grab_results(self):
        # grab initial conditions
        icfn = "%s_initial_conditions.h5" %(self.network_name)
        icf = h5py.File(icfn)
        self.icspecies = icf.keys()
        self.icnspecies = len(self.icspecies)
        self.icncells = icf[self.icspecies[0]].shape[0]
        self.icdata = dict(
            (s, np.empty((self.icncells), dtype="float64"))
            for s in self.icspecies)
        for s in self.icspecies:
            self.icdata[s][:] = icf[s][:]
        icf.close()

        # grab intermediate data
        fns = glob.glob("%s_intermediate_*.h5" % (self.network_name))
        if len(fns) == 0: raise RuntimeError
        f = h5py.File(fns[0])
        self.species = f.keys()
        print self.species
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
            self.t[i] = f["/"].attrs["time"]
            self.dt[i] = f["/"].attrs["timestep"]
            f.close()
        print self.t
        print "Loaded %s files" % (len(fns))

        # grab the final solution
                # grab initial conditions
        fsfn = "%s_solution.h5" %(self.network_name)
        fsf = h5py.File(fsfn)
        self.fsspecies = fsf.keys()
        self.fsnspecies = len(self.fsspecies)
        self.fsncells = fsf[self.fsspecies[0]].shape[0]
        self.fsdata = dict(
            (s, np.empty((self.fsncells), dtype="float64"))
            for s in self.fsspecies)
        for s in self.fsspecies:
            self.fsdata[s][:] = fsf[s][:]
        fsf.close()

    def plot_all(self):
        mpl.rcParams['axes.color_cycle'] = [list(clr) for clr in mpl.cm.jet(np.linspace(0,1,(self.nspecies - 1)))]
        mpl.rc('xtick', labelsize=14)
        mpl.rc('ytick', labelsize=14)
        mpl.rc('axes', labelsize=18)
        # plot the time evolution of a species for a given T
        plt.clf()
        fig = plt.figure(1, dpi=100)
        fig.set_size_inches(12,8)
        for j in range(self.ncells):
            i = 0
            for s in sorted(self.species):
                if s != 'ge':
                    i += 1
                    if s == 'T':
                        plt.loglog(self.t, self.data[s][j,:], '--', label='$T (K)$', lw=1.5)
                    else:
                        if s == 'de':
                            plt.loglog(self.t, self.data[s][j,:], '-.', label='$e^{-}$', lw=1.5, marker='o', markersize=5.0)
                        else:                        
                            plt.loglog(self.t, self.data[s][j,:], '-', label=s, lw=1.5)
            plt.xlabel("Time (s)")
            plt.ylabel("n$_{s}$")
            plt.ylim(1e-8, 2e7)
            #plt.xlim(1e4,1e12)
            ax = plt.axes()
            box = ax.get_position()
            ax.set_position([box.x0 - 0.05, box.y0, box.width * 0.9, box.height * 1.1])
            leg = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), ncol=2)
            #leg = plt.legend(loc = 'best', ncol=2)
            # leg.legendPatch.set_alpha(0.8)
            plt.savefig("%s_time_%04i.png" % (self.network_name, j))
            plt.clf()
        mpl.rcParams['axes.color_cycle'] = [list(clr) for clr in mpl.cm.spectral(np.linspace(0,1,(self.nspecies)))]
        plt.clf()
        for s in sorted(self.species):
            #if s != 'ge':
            if s == 'HII':
                plt.loglog(self.t, self.data[s][-1,:]/self.data[s][0,0], '-', label=s, lw=1.5, marker='x')
            else:
                plt.loglog(self.t, self.data[s][-1,:]/self.data[s][0,0], '-', label=s, lw=1.5)
        plt.xlabel("Time")
        plt.ylim(1e-4, 1e1)
        leg = plt.legend(loc = "best")
        leg.legendPatch.set_alpha(0.1)
        plt.savefig("%s_time_relative.png" % (self.network_name))
        
        # plot the initial conditions and the final solution
        mpl.rcParams['axes.color_cycle'] = [list(clr) for clr in mpl.cm.spectral(np.linspace(0,1,(self.nspecies - 2)))]
        plt.clf()
        for s in self.species:
            if s not in ('ge', 'T'):
                if s != 'de':
                    plt.loglog(self.icdata['ge'][:], self.icdata[s][:]/16.0, label=s, lw=1.5)
                else:
                    plt.loglog(self.icdata['ge'][:], self.icdata[s][:], label=s, lw=1.5)
        for s in self.species:
            if s not in ('ge', 'T'):
                plt.loglog(self.fsdata['ge'][:], self.fsdata[s][:], '--', lw=1.5)
        plt.xlabel("Temperature")
        plt.ylim(1e-3, 10)
        leg = plt.legend(loc = "best")
        leg.set_alpha(0.1)
        plt.savefig("%s_ic_final.png" %(self.network_name))
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("network_name", nargs=1)
    args = parser.parse_args()
    rp = ResultsPlotter(args.network_name[0])
    rp.plot_all()
