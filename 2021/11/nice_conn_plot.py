# Little script to just make a nice plot showing the connection lengths from MAFOT.
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import LimWallToolkit
import pickle
import matplotlib.path        as mplpath
import matplotlib as mpl


mafot_file1 = "/Users/zamperini/Documents/d3d_work/184527/mafot/lam_hires_tor0_conn+1.dat"
mafot_file2 = "/Users/zamperini/Documents/d3d_work/184527/mafot/lam_hires_tor0_conn-1.dat"
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/184527/184527_3500.pickle"
wall_path = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"
tor_angle = 0

# Copied from LimWallToolkit.
# Useful to check some of the data.
debug_dict = {}

# Load data into DataFrame.
print("Loading MAFOT runs...")
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
try:
    df1 = pd.read_csv(mafot_file1, skiprows=52, names=columns,
      delimiter="\t")
except FileNotFoundError:
    print("Error: Unable to find file: {}".format(mafot_file1))
    print("Exiting")
    sys.exit()
try:
    df2 = pd.read_csv(mafot_file2, skiprows=52, names=columns,
      delimiter="\t")
except FileNotFoundError:
    print("Error: Unable to find file: {}".format(mafot_file2))
    print("Exiting")
    sys.exit()

# Also read the file to pull out the number of R and Z coords.
with open(mafot_file1) as f:
    for line in f:
        if line[:10] == "# R-grid: ":
            numrs1 = int(line.split(":")[1])
        if line[:10] == "# Z-grid: ":
            numzs1 = int(line.split(":")[1])
            break
with open(mafot_file2) as f:
    for line in f:
        if line[:10] == "# R-grid: ":
            numrs2 = int(line.split(":")[1])
        if line[:10] == "# Z-grid: ":
            numzs2 = int(line.split(":")[1])
            break

# Reshape into 2D arrays. Reasonable assumption that both use the
# same R, Z.
if df1.shape != df2.shape:
    print("Error: Please use the same number of R and Z coordinates" + \
      " for both MAFOT runs.")
    #numrs1 = len(df1["R (m)"].unique())
    #numrs2 = len(df2["R (m)"].unique())
    #numzs1 = len(df1["Z (m)"].unique())
    #numzs2 = len(df2["Z (m)"].unique())
    print("            | # R's | # Z's |")
    print("mafot_file1 | {:5} | {:5} |".format(numrs1, numzs1))
    print("mafot_file2 | {:5} | {:5} |".format(numrs2, numzs2))
    print("Exiting")
    sys.exit()
r = df1["R (m)"].unique()
z = df1["Z (m)"].unique()
#r = df1["R (m)"][:numrs1]
#z = df1["Z (m)"].unique()
R, Z = np.meshgrid(r, z)

# Pull out each's connection length and pitch angles.
l1 = df1["Lconn (km)"].values * 1000  # km to m
l2 = df2["Lconn (km)"].values * 1000
p1 = df1["pitch angle"].values
p2 = df2["pitch angle"].values
L1 = l1.reshape(len(r), len(z))
L2 = l2.reshape(len(r), len(z))
pitch1 = p1.reshape(len(r), len(z))
pitch2 = p2.reshape(len(r), len(z))

# G-file needed for the equilibrium related info. This is not the actual
# gfile, but rather a pickled version from the OMFIT EFIT module. See
# comment at bottom of lwt_control_file.py file for instructions.
with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]
debug_dict["gfile"] = gfile

lwt = LimWallToolkit.LimWallToolkit()
wall = lwt.read_3d_wall(wall_path)
wall_coords = wall[tor_angle]

bbpath = mplpath.Path(list(zip(wall_coords[0], wall_coords[1])))
bbpath_mask = ~bbpath.contains_points(np.array(list(zip(gR.flatten(),
  gZ.flatten()))))
psin_masked = np.ma.masked_array(psin.flatten(),
  mask=bbpath_mask).reshape(psin.shape)
bbpath_mask2 = ~bbpath.contains_points(np.array(list(zip(R.flatten(),
  Z.flatten()))))
conn_masked = np.ma.masked_array((L1+L2).flatten(),
  mask=bbpath_mask2).reshape(L1.shape)

# Now actual plotting.
fig, ax = plt.subplots()

ax.set_aspect("equal")
ax.plot(wall_coords[0], wall_coords[1], color="k", lw=2)
ax.contour(gR, gZ, psin_masked, levels=[1], colors="k")
norm = mpl.colors.Normalize(vmin=0.05, vmax=50)
cont = ax.contourf(R, Z, conn_masked,
  cmap="inferno", vmin=0.05, vmax=50, levels=np.linspace(0.05, 50, 10), extend="max")
cbar = fig.colorbar(cont)
cbar.set_label("Connection Length (m)", fontsize=14)

fig.show()
