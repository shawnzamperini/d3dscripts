import oedge_plots
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

divimp_nc_path = "/Users/zamperini/Documents/d3d_work/167196/d3d-167196-mrc-shifted-nodrift.nc"
op = oedge_plots.OedgePlots(divimp_nc_path)

fig = op.plot_contour_polygon(dataname="DDLIMS", cmap="inferno", charge="all", vmin=1e13,
  vmax=3e14, scaling=op.absfac, normtype="log", cbar_label="W Density (m-3)")

ax = fig.axes[0]
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)

cbar_ticks = np.append(np.linspace(1e13, 1e14, 10), np.linspace(1e14, 3e14, 3))
#cbar_ticklabels = ["{:.0e}".format(l) for l in cbar_ticks]
cbar_ticklabels = []
for num in cbar_ticks:
    num_str = "{:.0e}".format(num)
    if num in [1e13, 1e14]:
        cbar_ticklabels.append(num_str)
    else:
        #cbar_ticklabels.append(num_str[0])
        cbar_ticklabels.append("")
cbar = fig.axes[1]
cbar = ax.collections[0].colorbar
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels(cbar_ticklabels)
