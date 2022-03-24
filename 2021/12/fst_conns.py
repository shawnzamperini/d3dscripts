# Script to plot the connection lengths for each CP side.
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.patches import Rectangle, Polygon


# All copied from check_lim.
ncpath = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_44.nc"
pol_idx = 21
nc = netCDF4.Dataset(ncpath)
cl = float(nc['CL'][:].data)
ca  = float(nc['CA'][:].data)
caw = float(nc['CAW'][:].data)
x = nc.variables['XOUTS'][:].data
y = nc.variables['YOUTS'][:].data
p = nc.variables["PS"][:].data
xkeep_min = np.nonzero(x)[0].min()
xkeep_max = np.nonzero(x)[0].max()
ykeep_min = np.nonzero(y)[0].min()
ykeep_max = np.nonzero(y)[0].max()
x = x[xkeep_min:xkeep_max]
y = y[ykeep_min:ykeep_max]
ykeep_cl = np.where(np.abs(y) < cl)[0]
y = y[ykeep_cl]
bounds1a = nc.variables["bounds_1a"][:].data
bounds2a = nc.variables["bounds_2a"][:].data
bounds1a_slice = bounds1a[pol_idx]
bounds2a_slice = bounds2a[pol_idx]

# Plotting.
fig, ax = plt.subplots(figsize=(5,4))

ax.step(bounds1a_slice, x, color="k", where="post")
ax.step(bounds2a_slice, x, color="k", where="post")

# Put a CP in the image.
width = 0.2
height = x.min() * -1
rect = Rectangle((0-width/2, -height), width, height, facecolor="k")
ax.add_patch(rect)

# Everything outside the volume grey.
poly1_coords = []; poly2_coords = []
for i in range(0, len(x)):
    if i < len(x)-1:
        poly1_coords.append([bounds1a_slice[i], x[i]])
        poly1_coords.append([bounds1a_slice[i+1], x[i]])
        poly2_coords.append([bounds2a_slice[i], x[i]])
        poly2_coords.append([bounds2a_slice[i+1], x[i]])
    else:
        poly1_coords.append([bounds1a_slice[i], x[i]])
        poly1_coords.append([bounds1a_slice[i]+10, x[i]])
        poly1_coords.append([bounds1a_slice[-1]+10, -0.2])
        poly1_coords.append([bounds1a_slice[0], -0.2])

        poly2_coords.append([bounds2a_slice[i], x[i]])
        poly2_coords.append([bounds2a_slice[i]-10, x[i]])
        poly2_coords.append([bounds2a_slice[-1]-10, -0.2])
        poly2_coords.append([bounds2a_slice[0], -0.2])
#poly1_coords = list(zip(bounds1a_slice, x))

poly1 = Polygon(poly1_coords, closed=True, facecolor="grey")
poly2 = Polygon(poly2_coords, closed=True, facecolor="grey")
ax.add_patch(poly1)
ax.add_patch(poly2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlabel("Parallel distance from probe (m)", fontsize=12)
ax.set_ylabel("Distance from probe tip (m)", fontsize=12)
ax.set_xlim([-10, 10])
ax.set_ylim([-0.11, 0.02])
ax.set_xticks(np.arange(-10, 15, 5))
yticks = np.arange(-0.10, 0.02, 0.02)
yticklabels = ["{:.2f}".format(-y) for y in yticks]
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

fig.tight_layout()
fig.show()
