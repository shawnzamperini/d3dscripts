import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# An example netCDF file. Doesn't matter too much which one.
ncpath = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_27b.nc"
nc = netCDF4.Dataset(ncpath)

# The value of the bounds gives its Z coordinate.
bounds1a = nc.variables["bounds_1a"][:].data.flatten()
bounds2a = nc.variables["bounds_2a"][:].data.flatten()

# Need X and P for each bounds location.
x = nc.variables['XOUTS'][:].data
xwids = nc.variables['XWIDS'][:].data
y = nc.variables['YOUTS'][:].data
p = nc.variables["PS"][:].data
pwids = nc.variables["PWIDS"][:].data
xkeep_min = np.nonzero(x)[0].min()
xkeep_max = np.nonzero(x)[0].max()
ykeep_min = np.nonzero(y)[0].min()
ykeep_max = np.nonzero(y)[0].max()
x = x[xkeep_min:xkeep_max]
xwids = xwids[xkeep_min:xkeep_max]
y = y[ykeep_min:ykeep_max]

X, P = np.meshgrid(x, p)
Xwids, Pwids = np.meshgrid(xwids, pwids)


# Width of bars.
dx = Xwids.flatten()
dy = Pwids.flatten()
z = np.zeros(len(bounds1a))


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.bar3d(X.flatten(), P.flatten(), z, dx, dy, bounds1a, shade=True)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

fig.show()
