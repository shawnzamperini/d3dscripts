# Script just a first-time glance at some of the results of running 3DLIM at
# various toroidal angles to look at how the impurity distibution changes.
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from LimWallToolkit import LimWallToolkit


# Some constants.
lim_path  = "/Users/zamperini/Documents/lim_runs/184527-tor30.nc"
r_origin3 = 2.295
z_origin3 = -0.188
wall_path3 = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"
along_coord3 = "R"
plot_surfaces = True
gfile_pickle_path3 = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"

# Load NC files and paths.
ncroot = "/users/zamperini/Documents/lim_runs/"
tors = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
ncs = []; ncpaths = []
for i in range(0, len(tors)):
    #ncpath = ncroot + "actual-184527-tor{}.nc".format(tors[i])
    ncpath = ncroot + "184527-tor{}.nc".format(tors[i])
    ncs.append(netCDF4.Dataset(ncpath))
    ncpaths.append(ncpath)

# Use the toolkit to get the goodies.
results = {}
lwt = LimWallToolkit()
for i in range(0, len(tors)):
    print("Angle {}".format(tors[i]))
    results[tors[i]] = lwt.plot_3dlim_on_wall(lim_path=ncpaths[i], tor_angle=tors[i],
      wall_path=wall_path3, plot_surfaces=True, gfile_pickle_path=gfile_pickle_path3,
      r_origin=r_origin3, z_origin=z_origin3, along_coord=along_coord3, show_plot=False)
    plt.close()

# All the 3DLIM results done on the same grid to simplify things. Calculate
# the difference in density from angle = 0.
diffs = {}
base_angle = 150
nz0 = results[base_angle]["masked_ddlim3"]
for i in range(0, len(tors)):
    nz = results[tors[i]]["masked_ddlim3"]
    diffs[tors[i]] = nz - nz0

    # Percent difference in nz from nz0.
    #diffs[tors[i]] = (nz - nz0) / nz0

machR = results[0]["lim_machRs"]
machZ = results[0]["lim_machZs"]

fig, axs = plt.subplots(3, 4, sharex=True, sharey=True)
for i in range(0, len(tors)):
    ax = axs.flatten()[i]
    lim = np.max(np.abs(diffs[tors[i]]))
    if lim == 0.0:
        levels = [-1, 0, 1]
    else:
        levels = np.linspace(-10, 10, 11)
    ax.pcolormesh(machR, machZ, diffs[tors[i]] + nz0, shading="auto", vmin=nz0.min(), vmax=nz0.max())
    #ax.contourf(machR, machZ, diffs[tors[i]], levels=levels, cmap="coolwarm", extend="both")
    wall_r = results[tors[i]]["wall_coords"][0]
    wall_z = results[tors[i]]["wall_coords"][1]
    ax.plot(wall_r, wall_z, color="k")
    ax.set_title(tors[i])
    ax.set_aspect("equal")
fig.show()

fig, axs = plt.subplots(3, 4, sharex=True, sharey=True)
for i in range(0, len(tors)):
    ax = axs.flatten()[i]
    R, Z = np.meshgrid(results[tors[i]]["lim_rbins"], results[tors[i]]["lim_pbins"])
    #ax.pcolormesh(R, Z, diffs[tors[i]], shading="auto", cmap="coolwarm",
    #  vmin=-0.00001, vmax=0.00001)
    lim = np.max(np.abs(diffs[tors[i]]))
    if lim == 0.0:
        levels = [-1, 0, 1]
    else:
        levels = np.linspace(-10, 10, 11)
    cont = ax.contourf(R, Z, diffs[tors[i]].T, cmap="coolwarm",
      levels=levels, extend="both")
    ax.set_title(tors[i])
    ax.set_xlim(R.max(), R.min())

fig.tight_layout()
cbar = fig.colorbar(cont, ax=axs)
cbar.set_label("C Density Difference (arbitrary)")
fig.show()

# A More concise plot of density at a point plotted as function of toroidal angle.
test_R = 2.34
test_Z = 0.00
test_dist = np.full(machR.shape, 9999.0)
for i in range(0, machR.shape[0]):
    for j in range(0, machR.shape[1]):
        tmp = np.sqrt(np.square(machR[i, j] - test_R) + np.square(machZ[i, j] - test_Z))
        test_dist[i, j] = tmp
test_idx = np.where(test_dist == test_dist.min())
single_diffs = []
for i in range(0, len(tors)):
    single_diffs.append(diffs[tors[i]][test_idx].data[0])
single_diffs = np.array(single_diffs)
fig, ax = plt.subplots()
ax.plot(tors, single_diffs * 100, marker="o", ms=5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_title("% Difference from nC @ {}".format(base_angle), fontsize=16)
ax.set_xlabel("Toroidal Angle", fontsize=14)
#ax.set_ylim([25, 150])
ax.axhline(0.0, color="k")
fig.tight_layout()
fig.show()
