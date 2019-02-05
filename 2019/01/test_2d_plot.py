import netCDF4
import matplotlib.pyplot  as plt
import numpy              as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches


# Constants
iz_state = 5
r_min = -0.01
r_max = 0.0
cp_width = 3.0

# Load netcdf file.
netcdf = netCDF4.Dataset('colprobe-test-m2.nc')

# Get positions of the center of each bin.
xs       = netcdf.variables['XS'][:].data
xwids    = netcdf.variables['XWIDS'][:].data
rad_locs = xs-xwids/2.0
ps       = netcdf.variables['PS'][:].data
pwids    = netcdf.variables['PWIDS'][:].data
pol_locs = ps-pwids/2.0

# Same for parallel bins, but get rid of the extra zeros.
ys       = netcdf.variables['YS'][:].data
ywids    = netcdf.variables['YWIDS'][:].data
par_locs = ys-ywids/2.0

# Also mirror the par_locs to cover both sides (i.e. from -L to L instead of 0 to L).
# Need to add a zero in the as the middle point, hence two appends.
par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

# Load ddlim3 variable array.
ddlim3 = netcdf.variables['DDLIM3'][:, iz_state, :, :].data

sum_range = np.where(np.logical_and(rad_locs>r_min, rad_locs<r_max))[0]
summed_ddlim3 = ddlim3[:,:,sum_range].sum(axis=2)

X, Y = np.meshgrid(par_locs, pol_locs)
Z = summed_ddlim3
fig = plt.figure()
ax = fig.add_subplot(111)
cont = ax.contour(X, Y, Z)
cbar = fig.colorbar(cont)
ax.set_xlim([-10, 10])
cp = patches.Rectangle((-0.2,-0.015), width=0.4, height=0.03, color='k')
ax.add_patch(cp)
fig.show()
