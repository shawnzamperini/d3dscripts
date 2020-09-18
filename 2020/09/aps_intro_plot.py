import oedge_plots
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


# Some constants.
ring_end = 50
pfz_start = 72
alpha = 0.5

# Load case.
ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-bkg-011f.nc'
op = oedge_plots.OedgePlots(ncpath)

# Better wall data.
keep_idx = np.where(np.logical_and(op.rvesm[0] != 0, op.zvesm[0] != 0))
rvesm = np.append(op.rvesm[0][keep_idx], op.rvesm[0][keep_idx][0])
zvesm = np.append(op.zvesm[0][keep_idx], op.zvesm[0][keep_idx][0])
rvesm[62:71] = np.full(len(rvesm[62:71]), rvesm[62])
zvesm[62:71] = np.full(len(zvesm[62:71]), zvesm[62])

# Load rings data, assign values according to color wanted.
rings = op.read_data_2d("KTEBS", scaling="Ring")
near_sol = np.logical_and(rings >= ring_end, rings < pfz_start)
far_sol = np.logical_and(rings >= 19, rings < ring_end)
core = rings < 19
pfz = rings > pfz_start
rings[near_sol] = 0
rings[pfz] = 0.5
rings[core] = 0.5
rings[far_sol] = 1

# Custom colormap.
cmap = ListedColormap([(44/255, 160/255, 44/255, alpha), (0, 0, 0, 0), (214/255, 39/255, 40/255, alpha)])
#cmap = cm.get_cmap('coolwarm', 3)

# Plot.
op.plot_contour_polygon("KTEBS", cmap=cmap, own_data=rings, vmin=0, vmax=1,
                        wall_data=(rvesm, zvesm), no_core=True)
