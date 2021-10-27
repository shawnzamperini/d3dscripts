# Script to make a plot of the 2D Ti and C distributions.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge")
import oedge_plots


ncpath = "/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-007.nc"
op = oedge_plots.OedgePlots(ncpath)

def plot(tag, vmin, vmax):

    if tag == "DDLIMS":
        charge     = "all"
        scaling    = op.absfac
        cbar_label = "C13 Density (m-3)"
    else:
        charge     = None
        scaling    = None
        cbar_label = "Ti (eV)"

    fig = op.plot_contour_polygon(tag, charge=charge, cmap="magma", vmin=vmin,
      vmax=vmax, normtype="log", scaling=scaling, cbar_label=cbar_label)

    # Remove spines.
    fig.axes[0].spines["top"].set_visible(False)
    fig.axes[0].spines["right"].set_visible(False)
    fig.axes[0].spines["bottom"].set_visible(False)
    fig.axes[0].spines["left"].set_visible(False)
    fig.axes[0].set_xticks([])
    fig.axes[0].set_yticks([])
    fig.axes[0].set_xlabel("")
    fig.axes[0].set_ylabel("")

    fig.show()

plot("DDLIMS", vmin=1e14, vmax=1e16)
