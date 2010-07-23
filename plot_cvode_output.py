import matplotlib;matplotlib.use("Agg");import pylab
import h5py
import numpy as na

YINS = 365*24*3600

def read_datasets(f):
    print "Reading from", f
    vals = []
    times = []
    rho = f["/"].attrs["rho"][:]
    for dsname in sorted(f["/"]):
        times.append(f[dsname].attrs["t"])
        vals.append(f[dsname][:])
    return rho, na.array(vals), na.array(times)

def get_varnames(f):
    print "Getting varnames from", f
    id_lookup = dict([(int(a), "".join(v.view("c")))
                    for a, v in f["/"].attrs.items()
                    if a not in ["rho"] ])
    var_lookup = dict( [(b,a) for a, b in id_lookup.items()] )
    return id_lookup, var_lookup

def next_ls():
    while 1:
        for ls in ["-", "-.", ":", "--"]:
            yield ls
lss = next_ls()

def compute_temperature(vals, rho):
    rho_cgs = rho[0] * 1.67e-24
    p2d = (5./3. - 1.0) * rho_cgs * vals["ge"]
    tgas = 0.0
    tgas += vals["HI"] + vals["HII"] + vals["HM"] + vals["de"]
    tgas += (vals["H2I"] + vals["H2II"])/2.0
    tgas += (vals["HeI"] + vals["HeII"] + vals["HeIII"])/4.0
    return p2d / (1.380e-16 * tgas)

def plot_vars(vals, rho, times, var_lookup, id_lookup):
    lvals = dict( [(var, vals[:,0,vid]) for var, vid in sorted(var_lookup.items())])
    pylab.clf()
    fig = pylab.gcf()
    ax = pylab.subplot(2,1,1)
    h = {}
    for var, vid in sorted(var_lookup.items()):
        if var == "ge": continue
        h[var] = ax.loglog(times/YINS, vals[:,0,vid], lw = 1.5, ls = lss.next())
    ax.set_ylim( rho.max()*1e-16, rho.max()*2.0 )
    labels, handles = zip(*h.items())
    pylab.figlegend(handles, labels, "upper right",
                    prop = dict(size = "x-small"))
    #ax.legend(loc = "upper left")
    ax.set_ylabel("amu / cc")
    ax.xaxis.set_ticklabels(())
    ax = pylab.subplot(2,1,2)
    T = compute_temperature(lvals, rho)
    ax.semilogx(times/YINS, T)
    ax.set_ylabel("T")
    ax.set_xlabel("years")
    pylab.savefig("output.pdf")

if __name__ == "__main__":
    f = h5py.File("primordial_output.h5")
    id_lookup, var_lookup = get_varnames(f)
    rho, vals, times = read_datasets(f)
    plot_vars(vals, rho, times, var_lookup, id_lookup)
