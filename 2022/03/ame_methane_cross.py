# Script to plot the cross section of a methane shot using just the EFIT data.
import pickle
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, Point
from tqdm import tqdm


# Load gfile, pull out relevant arrays.
pickle_path = "/Users/zamperini/Documents/d3d_work/lwt/184527/184527_3000.pickle"
with open(pickle_path, "rb") as f:
    gfile = pickle.load(f)

wall_r = gfile["RLIM"]
wall_z = gfile["ZLIM"]
raxis = gfile["RMAXIS"]
zaxis = gfile["ZMAXIS"]
r = gfile["R"]
z = gfile["Z"]
psin = gfile["PSIRZ_NORM"]
R, Z = np.meshgrid(r, z)

# Make a Polygon out of the wall. Identify points which exist in it.
wall = Polygon(zip(wall_r, wall_z))
inside = np.full(R.shape, False)
for i in tqdm(range(0, R.shape[0])):
    for j in range(0, Z.shape[1]):
        p = Point(R[i, j], Z[i, j])
        if wall.contains(p):
            inside[i, j] = True
psin = np.ma.masked_where(~inside, psin)

fig, ax = plt.subplots(figsize=(4, 5))
ax.plot(wall_r, wall_z, color="k", lw=2)
ax.set_aspect("equal")
ax.scatter(raxis, zaxis, marker="+", s=60, color="k")
ax.contour(R, Z, psin, colors="k", levels=[1], linewidths=2)
ax.contour(R, Z, psin, colors="k", levels=np.linspace(1.0, 1.2, 6))
ax.axis("off")
fig.tight_layout()
fig.show()
