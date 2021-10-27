# Script just a first-time glance at some of the results of running 3DLIM at
# various toroidal angles to look at how the impurity distibution changes.
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from LimWallToolkit import LimWallToolkit


# Some constants.
rmrsep_origin = 0.037
p_to_z = [0.9766, -0.188]
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"
wall_path = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"

# Load NC files and paths.
ncroot = "/users/zamperini/Documents/lim_runs/"
tors = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
ncs = []; ncpaths = []
for i in range(0, len(tors)):
    ncpath = ncroot + "actual-184527-tor{}.nc".format(tors[i])
    ncs.append(netCDF4.Dataset(ncpath))
    ncpaths.append(ncpath)

# Use the toolkit to get the goodies.
results = {}
lwt = LimWallToolkit()
for i in range(0, len(tors)):
    print("Angle {}".format(tors[i]))
    results[tors[i]] = lwt.plot_3dlim_on_wall(lim_path=ncpaths[i], tor_angle=tors[i],
      rmrsep_origin=rmrsep_origin, p_to_z=p_to_z, wall_path=wall_path,
      plot_surfaces=True, gfile_pickle_path=gfile_pickle_path, show_plot=False)

# All the 3DLIM results done on the same grid to simplify things. Calculate
# the difference in density from angle = 0.
diffs = {}
nz0 = results[0]["imp_density"]
for i in range(0, len(tors)):
    nz = results[tors[i]]["imp_density"]
    diffs[tors[i]] = nz - nz0

R = results[0]["R_lim_2d"]
Z = results[0]["zs_2d"]

fig, axs = plt.subplots(3, 4, sharex=True, sharey=True)
for i in range(0, len(tors)):
    ax = axs.flatten()[i]
    ax.pcolormesh(R, Z, diffs[tors[i]], shading="auto")
    wall_r = results[tors[i]]["wall"][0]
    wall_z = results[tors[i]]["wall"][1]
    ax.plot(wall_r, wall_z, color="k")
    ax.set_title(tors[i])
fig.show()

fig, axs = plt.subplots(3, 4, sharex=True, sharey=True)
for i in range(0, len(tors)):
    ax = axs.flatten()[i]
    R, Z = np.meshgrid(results[tors[i]]["rad_locs"], results[tors[i]]["pol_locs"])
    #ax.pcolormesh(R, Z, diffs[tors[i]], shading="auto", cmap="coolwarm",
    #  vmin=-0.00001, vmax=0.00001)
    cont = ax.contourf(R, Z, diffs[tors[i]], cmap="coolwarm",
      levels=np.linspace(-0.00005, 0.00005, 10), extend="both")
    ax.set_title(tors[i])
    ax.set_xlim(R.max(), R.min())

fig.tight_layout()
cbar = fig.colorbar(cont, ax=axs)
cbar.set_label("C Density Difference (arbitrary)")
fig.show()
