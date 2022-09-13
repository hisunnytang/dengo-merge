import matplotlib.pyplot as plt
import numpy as np
import yt
from matplotlib.colors import LogNorm, SymLogNorm
from yt.visualization.api import get_multi_plot


def plot_multi_sliceplot(
    dengo_ds,
    grackle_ds,
    fields,
    width,
    rule_scale,
    normal="z",
    cmaps=["bds_highcontrast", "plasma", "cmyt.pastel", "cmyt.dusk"],
    titles=None,
    dengo_center="m",
    grackle_center="m",
    fsize=512,
    use_symlog=[False, False, False, False],
    scale_color="white",
):
    orient = "horizontal"
    nfields = len(fields)

    fig, axes, colorbars = get_multi_plot(nfields, 2, colorbar=orient, bw=4)

    slc_dengo = yt.SlicePlot(
        dengo_ds,
        normal,
        center=dengo_center,
        fields=fields,
    )
    slc_grackle = yt.SlicePlot(
        grackle_ds,
        normal,
        center=grackle_center,
        fields=fields,
    )

    slc_dengo_frb = slc_dengo.data_source.to_frb(width, fsize)
    slc_grackle_frb = slc_grackle.data_source.to_frb(width, fsize)

    plots = []
    # make dengo plots
    i = 0
    for ax_dengo, ax_grackle, f, cmap in zip(axes[0], axes[1], fields, cmaps):
        ax_dengo.xaxis.set_visible(False)
        ax_dengo.yaxis.set_visible(False)

        ax_grackle.xaxis.set_visible(False)
        ax_grackle.yaxis.set_visible(False)

        slc = np.array(slc_dengo_frb[f])
        slc_min, slc_max = slc.min(), slc.max()

        if ((slc_min > 0) * (slc_max > 0)) and not use_symlog[i]:
            bnd = [slc_min, slc_max]
            norm = LogNorm()
        else:
            abs_max = np.max([np.abs(slc_min), np.abs(slc_max)])
            abs_min = np.min([np.abs(slc_min), np.abs(slc_max)])
            bnd = [-abs_max, abs_max]

            norm = SymLogNorm(linthresh=abs_max * 1e-2)

        im_dengo = ax_dengo.imshow(slc, origin="lower", norm=norm)

        slc = np.array(slc_grackle_frb[f])
        im_grackle = ax_grackle.imshow(slc, origin="lower", norm=norm)

        im_dengo.set_clim(bnd)
        im_dengo.set_cmap(cmap)

        im_grackle.set_clim(bnd)
        im_grackle.set_cmap(cmap)

        plots.append(im_grackle)

        i += 1

    # Annotate
    pixel_size = (rule_scale / width).v * fsize

    axes[1][-1].plot(
        [0.8 * fsize, 0.8 * fsize + pixel_size],
        [0.1 * fsize, 0.1 * fsize],
        color=scale_color,
    )
    axes[1][-1].annotate(
        f"{int(rule_scale.v.item())} {rule_scale.units}",
        xy=(0.8 * fsize, 0.05 * fsize),
        color=scale_color,
        fontsize=10,
    )

    if titles is None:
        titles = fields
    for p, cax, t in zip(plots, colorbars, titles):
        print(t)
        cbar = fig.colorbar(p, cax=cax, orientation=orient)
        cbar.set_label(t)

    return fig
