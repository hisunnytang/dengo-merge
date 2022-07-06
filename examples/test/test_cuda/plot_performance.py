import h5py
import matplotlib.pyplot as plt
import numpy as np


def plot_scaling():
    cuda_scaling = h5py.File("cuda_solver/scaling.h5", "r")
    suitesparse_scaling = h5py.File("suitesparse/scaling.h5", "r")

    cuda_times = cuda_scaling["time"][:]
    ss_times = suitesparse_scaling["time"][:]
    cells = cuda_scaling["ncells"][:]

    print(f"cuda:        {cuda_times/cells}")
    print(f"suitesparse: {ss_times/cells}")

    f, ax = plt.subplots()
    ax.loglog(cells, cuda_times / cells, label="Cuda")
    ax.loglog(cells, ss_times / cells, label="SuiteSparse-OpenMP")

    ax2 = ax.twinx()
    ax2.semilogx(cells, ss_times / cuda_times, color="r", ls="--", label="ratio")

    ax.set_xlabel("Number of Cells")
    ax.set_ylabel("Runtime/ Cell (s)")
    ax2.set_ylabel("Runtime Ratio")

    # ax.legend(loc=0)
    # ax2.legend(loc=2)
    ax.legend(bbox_to_anchor=(0.6, 0.1), loc="lower left", frameon=False)
    ax2.legend(bbox_to_anchor=(0.6, 0.05), loc="lower left", frameon=False)
    plt.tight_layout()
    f.savefig("comparison.png")

    cuda_scaling.close()
    suitesparse_scaling.close()


def plot_performance():

    cuda_perf = h5py.File("cuda_solver/performance.h5", "r")
    suitesparse_perf = h5py.File("suitesparse/performance.h5", "r")

    cuda_times = cuda_perf["time"][:]
    ss_times = suitesparse_perf["time"][:]

    d = cuda_perf["density"][:]
    t = cuda_perf["temperature"][:]

    d2d = d.reshape(9, 9)
    t2d = t.reshape(9, 9)

    print(f"cuda:        {cuda_times}")
    print(f"suitesparse: {ss_times}")

    f, ax = plt.subplots(1, 1)
    ratio = np.log10(ss_times / cuda_times)
    ax.hist(ratio, label="ratio")
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect((xright - xleft) / (ytop - ybottom))
    ax.set_xlabel("$\log \quad t_{ss}/ t_{cuda}$")
    ax.set_ylabel("Count")
    plt.tight_layout()
    f.savefig("grid_comparison_hist.png")

    plt.clf()

    f, ax = plt.subplots(1, 1)
    logratio2d = ratio.reshape(9, 9)
    c = ax.pcolormesh(
        d2d, t2d, logratio2d, cmap="RdBu", vmin=min(ratio), vmax=max(ratio)
    )
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Density ($cm^{-3}$)")
    ax.set_ylabel("Temperature ($K$)")

    xleft, xright = np.log(ax.get_xlim())
    ybottom, ytop = np.log(ax.get_ylim())
    ax.set_aspect((xright - xleft) / (ytop - ybottom))

    f.colorbar(c, ax=ax, label="$\log \quad t_{ss}/ t_{cuda}$")

    plt.tight_layout()
    f.savefig("grid_comparison.png")

    cuda_perf.close()
    suitesparse_perf.close()


if __name__ == "__main__":
    plot_performance()
    plot_scaling()
