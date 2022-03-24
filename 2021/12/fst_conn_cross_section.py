# Script to mat a set of plots of the connection length contours.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from matplotlib import ticker, cm
from matplotlib.colors import LogNorm
from LimWallToolkit import LimWallToolkit



# Load in MAFOT data for 167196.
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
mafot_file1 = "/Users/zamperini/Documents/d3d_work/167196/lam_hires_tor240_conn+1.dat"  # CCW, Reverse =
mafot_file2 = "/Users/zamperini/Documents/d3d_work/167196/lam_hires_tor240_conn-1.dat"  # CW,  Reverse =
df1 = pd.read_csv(mafot_file1, skiprows=52, names=columns, delimiter="\t")
df2 = pd.read_csv(mafot_file2, skiprows=52, names=columns, delimiter="\t")
r = df1["R (m)"].unique()
z = df1["Z (m)"].unique()
R, Z = np.meshgrid(r, z)

# Pull out each's connection length.
l1 = df1["Lconn (km)"].values * 1000  # km to m
l2 = df2["Lconn (km)"].values * 1000
L1 = l1.reshape(len(r), len(z))
L2 = l2.reshape(len(r), len(z))
L = L1 + L2

# Masked length data to avoid core data.
L = np.ma.masked_where(np.logical_and(L>2000, L<=0), L)

# Load in gfile for the separatrix.
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/167196/167196_3500.pickle"
with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]

# Load wall coordinates.
wall_path = "/Users/zamperini/Documents/d3d_work/lwt/167196/mafot_3D_wall.dat"
lwt = LimWallToolkit()
wall = lwt.read_3d_wall(wall_path)
wall_coords = wall[240]

# Get the lengths of field lines at a value close to MiMES Z location.
mimes_z = -0.188
mimes_zidx = np.where(np.abs(Z-mimes_z) == np.abs((Z-mimes_z)).min())
mimes_L1 = L1[mimes_zidx]
mimes_L2 = L2[mimes_zidx]
mimes_R  = R[mimes_zidx]

# Plotting commands.
levels = np.geomspace(0.1, 100, 10)
fig, ax = plt.subplots()

ax.plot(wall_coords[0], wall_coords[1], color="k", zorder=3)
cont = ax.contourf(R, Z, L, cmap="inferno", norm=LogNorm(), vmin=0.1, levels=levels)
cbar = fig.colorbar(cont)
cbar.ax.minorticks_off()
#cbar.set_labels()

# Plot separatrix.
ax.contour(gR, gZ, psin, levels=[1], colors="k", linewidths=3, zorder=40)

cbar.set_label("Length of Field Line (m)", fontsize=12)
ax.set_aspect("equal")

fig.tight_layout()
fig.show()

# Another plot of MiMES Lconn.
fig, ax = plt.subplots(figsize=(5,4))
ax.plot(mimes_R, mimes_L1, lw=3, label="Outer Target Facing", color="tab:purple")
ax.plot(mimes_R, mimes_L2, lw=3, label="Inner Target Facing", color="tab:red")
ax.set_yscale("log")
ax.legend(fontsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim([2.20, 2.40])
ax.set_xticks(np.arange(2.20, 2.45, 0.05))
ax.set_ylim([0.01, 1000])
ax.set_xlabel("R (m)", fontsize=12)
ax.set_ylabel("Length of Field Line (m)", fontsize=12)
fig.tight_layout()
fig.show()
