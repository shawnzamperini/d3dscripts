import oedge_plots
import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams["font.family"] = "Century Gothic"
#plt.rc('axes', unicode_minus=False)

unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}.nc".format(167277, 167277, "006")
op = oedge_plots.OedgePlots(unf_path)
fig = op.plot_contour_polygon("DDLIMS", charge="all", cmap="magma", vmin=0.001, vmax=0.03, normtype="log", cbar_label="W Density (arbitrary)")

#fig.axes[1].remove()
ax = fig.axes[0]
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="both", which="both", labelleft=False, labelbottom=False, left=False, bottom=False)
