# Plot the wall coordinates at a given angle with the equlibrium.
from LimWallToolkit import LimWallToolkit
import matplotlib.pyplot as plt
import pickle
import numpy as np


# First load in the wall coordinates.
wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
lwt = LimWallToolkit()
wall = lwt.read_3d_wall(wall_path)

# Then load in the gfile data.
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/lwt/187111/187111_3300.pickle"
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/167196/167196_3500.pickle"
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/930116_2110.pickle"
with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]


# Plotting.
fig, ax = plt.subplots()

ax.contour(gR, gZ, psin, levels=np.linspace(1.1, 1.4, 3), colors="k", linewidths=1)
ax.contour(gR, gZ, psin, levels=[1.0], linewidths=3, colors="k")

ax.plot(wall[242][0], wall[242][1], lw=2, color="k")
ax.set_aspect("equal")

fig.tight_layout()
fig.show()
