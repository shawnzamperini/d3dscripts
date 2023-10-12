# Script to use the EFIT to calculate Lpol into the windowed region.
import pickle
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


with open("/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle", "rb") as f:
    gfile = pickle.load(f)

R, Z = np.meshgrid(gfile["R"], gfile["Z"])

# Will use the contour plot to get the contour line we want.
fig, ax1 = plt.subplots()
cont = ax1.contour(R, Z, gfile["PSIRZ_NORM"], levels=[1.0, 1.13], colors="k")
ax1.plot(gfile["RLIM"], gfile["ZLIM"], color="k")
ax1.set_aspect("equal")

# We set it up so the second line is the contour we are interested in. Extract the coordinates, mask to just values in
# the vessel.
coords = cont.collections[1].get_paths()[1].vertices
wall = Polygon(zip(gfile["RLIM"], gfile["ZLIM"]))

# We can use Z = 1.1 m as a cutoff.
mask = np.array([wall.contains(Point(c)) and c[1] < 1.1 for c in coords])
coords_in = coords[mask]
ax1.scatter(coords_in[:, 0], coords_in[:, 1], color="r")

fig.tight_layout()
fig.show()

# At this point we have our line, all that remains is to calculate the length of it.
lpol = 0
for i in range(1, len(coords_in)):
    dr = coords_in[i, 0] - coords_in[i-1, 0]
    dz = coords_in[i, 1] - coords_in[i-1, 1]
    lpol += np.sqrt(np.square(dr) + np.square(dz))
print("Lpol = {:.2f} m".format(lpol))