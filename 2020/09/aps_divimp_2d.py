import matplotlib.pyplot as plt
import oedge_plots
import numpy as np


# Load run data.
#ncpath1 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-025a.nc'
#cmap = "Oranges"
cmap = "gnuplot"
#ncpath1 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-025b.nc'
ncpath1 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-026d.nc'
#cmap = "Blues"
#cmap = "PuBuGn"
op1 = oedge_plots.OedgePlots(ncpath1)

# Better wall data.
keep_idx = np.where(np.logical_and(op1.rvesm[0] != 0, op1.zvesm[0] != 0))
rvesm = np.append(op1.rvesm[0][keep_idx], op1.rvesm[0][keep_idx][0])
zvesm = np.append(op1.zvesm[0][keep_idx], op1.zvesm[0][keep_idx][0])
rvesm[62:71] = np.full(len(rvesm[62:71]), rvesm[62])
zvesm[62:71] = np.full(len(zvesm[62:71]), zvesm[62])

# 2D impurity distribution plots with some tweaking.
fig = op1.plot_contour_polygon("DDLIMS", charge="all", cmap=cmap, vmin=0.001,
                               vmax=0.075, wall_data=(rvesm, zvesm), lut=11)
ax = fig.axes[0]
cbar = fig.axes[1]
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticks([])
ax.set_yticks([])
cbar.set_ylabel("W Density (normalized)")
