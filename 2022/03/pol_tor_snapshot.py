# This script just takes the type of plots in conn_psin_gif.py and does a
# side by side comparison at a given psin.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from LimWallToolkit import LimWallToolkit
from matplotlib.colors import LogNorm, Normalize
import matplotlib.tri as tri


# Inputs
shot = 167196
plot_psin = 1.245

# Set the correct paths for each case.
pol_wall_path = "/Users/zamperini/Documents/d3d_work/lwt/167196/mafot_3d_wall.dat"
tor_wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
if shot == 167196:
    pol_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_with_pol_lims.dat"
    pol_ncoords_path = "ncoords_167196_pol.pickle"
    tor_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_with_tor_lims.dat"
    tor_ncoords_path = "ncoords_167196_tor.pickle"
    gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/167196/167196_3500.pickle"

# Load in the mafot data into a dataframe.
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
pol_df = pd.read_csv(pol_mafot_path, skiprows=52, names=columns, delimiter="\t")
tor_df = pd.read_csv(tor_mafot_path, skiprows=52, names=columns, delimiter="\t")
pol_r = pol_df["R (m)"].values
pol_z = pol_df["Z (m)"].values
pol_l = pol_df["Lconn (km)"].values * 1000
tor_r = tor_df["R (m)"].values
tor_z = tor_df["Z (m)"].values
tor_l = tor_df["Lconn (km)"].values * 1000

# Get the dictionary with the number of coordinates per flux tube.
with open(pol_ncoords_path, "rb") as f:
    pol_ncoords = pickle.load(f)
with open(tor_ncoords_path, "rb") as f:
    tor_ncoords = pickle.load(f)

# Load in the wall coordinates and gfile.
lwt = LimWallToolkit()
pol_wall = lwt.read_3d_wall(pol_wall_path)
tor_wall = lwt.read_3d_wall(tor_wall_path)

with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin_rz = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]

# Limit to just the right of the OMP.
raxis = gfile["RMAXIS"]
outboard = gfile["R"] > raxis
gR_out = gR[:, outboard]
gZ_out = gZ[:, outboard]
psin_rz_out = psin_rz[:, outboard]

# The coordinates are printed one psin at a time, [0-360] before moving onto
# the next psin.
n = 0
minz = 999; maxz = -999
minl = 999; maxl = -999
pol_conns = {}
for psin, ncoord in pol_ncoords.items():
    tmp_l = []; tmp_r = []; tmp_z = []; tmp_d = []
    for deg in range(0, 360, 2):
        for i in range(0, ncoord):
            tmp_l.append(pol_l[n])
            tmp_r.append(pol_r[n])
            tmp_z.append(pol_z[n])
            tmp_d.append(deg)

            # Keep track of minimum and maximum values for plots.
            if pol_z[n] < minz:
                minz = pol_z[n]
            if pol_z[n] > maxz:
                maxz = pol_z[n]
            if pol_l[n] < minl:
                minl = pol_l[n]
            if pol_l[n] > maxl:
                maxl = pol_l[n]

            n += 1
    pol_conns[psin] = {"l":np.array(tmp_l), "r":np.array(tmp_r),
        "z":np.array(tmp_z), "d":np.array(tmp_d)}

n = 0
tor_conns = {}
for psin, ncoord in tor_ncoords.items():
    tmp_l = []; tmp_r = []; tmp_z = []; tmp_d = []
    for deg in range(0, 360, 2):
        for i in range(0, ncoord):
            tmp_l.append(tor_l[n])
            tmp_r.append(tor_r[n])
            tmp_z.append(tor_z[n])
            tmp_d.append(deg)

            # Keep track of minimum and maximum values for plots.
            if tor_z[n] < minz:
                minz = tor_z[n]
            if tor_z[n] > maxz:
                maxz = tor_z[n]
            if tor_l[n] < minl:
                minl = tor_l[n]
            if tor_l[n] > maxl:
                maxl = tor_l[n]

            n += 1
    tor_conns[psin] = {"l":np.array(tmp_l), "r":np.array(tmp_r),
        "z":np.array(tmp_z), "d":np.array(tmp_d)}

# For the colorbar.
norm = Normalize(minl, maxl)
levels = np.arange(0, maxl, 2)

# Find nearest psin value in the dictionary.
psins = np.array(list(pol_conns.keys()))
near_psin = psins[np.argmin(np.abs(psins - plot_psin))]
print("Nearest psin = {:.4}".format(near_psin))

# Set zeros to something outside the range so they don't show up.
pol_conns[near_psin]["l"][pol_conns[near_psin]["l"]==0] = -999
tor_conns[near_psin]["l"][tor_conns[near_psin]["l"]==0] = -999

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1.5, 3]}, figsize=(9,6))

flux_bounds = [psins[0], psins[-1]]
ax1.contour(gR, gZ, psin_rz, levels=flux_bounds, colors="k", linewidths=1)
ax1.contour(gR, gZ, psin_rz, levels=[1.0], linewidths=2, colors="k")
ax1.contour(gR_out, gZ_out, psin_rz_out, levels=[near_psin], colors="tab:red")
ax1.plot(pol_wall[242][0], pol_wall[242][1], lw=2, color="k")
ax1.set_aspect("equal")
ax1.axis("off")

tri2 = ax2.tricontourf(pol_conns[near_psin]["d"], pol_conns[near_psin]["z"],
    pol_conns[near_psin]["l"], cmap="inferno", norm=norm, levels=levels,
    extend="max")
ax2.set_facecolor("grey")
ax2.set_ylim(minz, maxz)
ax2.set_ylabel("\n\nZ (m)", fontsize=12)
ax2.set_xticks([])
cbar = fig.colorbar(tri2, ax=[ax2, ax4], ticks=levels)
cbar.set_label("Connection Length (m)", fontsize=12)

ax3.contour(gR, gZ, psin_rz, levels=flux_bounds, colors="k", linewidths=1)
ax3.contour(gR, gZ, psin_rz, levels=[1.0], linewidths=2, colors="k")
ax3.contour(gR_out, gZ_out, psin_rz_out, levels=[near_psin], colors="tab:red")
ax3.plot(tor_wall[242][0], tor_wall[242][1], lw=2, color="k")
ax3.set_aspect("equal")
ax3.axis("off")

tri4 = ax4.tricontourf(tor_conns[near_psin]["d"], tor_conns[near_psin]["z"],
    tor_conns[near_psin]["l"], cmap="inferno", norm=norm, levels=levels,
    extend="max")
ax4.set_facecolor("grey")
ax4.set_ylim(minz, maxz)
ax4.set_xlabel("Toroidal Angle", fontsize=12)
ax4.set_ylabel("\n\nZ (m)", fontsize=12)

#fig.tight_layout()
fig.show()
