# Plot the wall coordinates at a given angle with the equlibrium.
from LimWallToolkit import LimWallToolkit
import matplotlib.pyplot as plt
import pickle
import numpy as np
from shapely.geometry import Polygon, Point


# First load in the wall coordinates.
wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
lwt = LimWallToolkit()
wall = lwt.read_3d_wall(wall_path)

# Then load in the gfile data.
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/lwt/187111/187111_3300.pickle"
#gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/167196/167196_3500.pickle"
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/930116_2110.pickle"
with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]


# Use most current wall.
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/lwt/187111/187111_3300.pickle"
with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
rwall = gfile["RLIM"]
zwall = gfile["ZLIM"]

# Mask points outside of the wall.
mask = np.full(gR.shape, True)
wall_poly = Polygon(list(zip(rwall, zwall)))
for i in range(gR.shape[0]):
    for j in range(gZ.shape[1]):
        p = Point(gR[i,j], gZ[i,j])
        if wall_poly.contains(p):
            mask[i,j] = False

gRm = np.ma.masked_where(mask, gR)
gZm = np.ma.masked_where(mask, gZ)
psinm = np.ma.masked_where(mask, psin)

# Plotting.
fig, ax = plt.subplots()

ax.contour(gRm, gZm, psinm, levels=np.geomspace(1.0, 1.2, 5), colors="k", linewidths=1)
ax.contour(gRm, gZm, psinm, levels=[1.0], linewidths=3, colors="k")

#ax.plot(wall[242][0], wall[242][1], lw=2, color="k")
ax.plot(rwall, zwall, lw=2, color="k")
ax.set_aspect("equal")

fig.tight_layout()
fig.show()
