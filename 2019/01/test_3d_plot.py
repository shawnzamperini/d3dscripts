import netCDF4
import matplotlib.pyplot  as plt
import numpy              as np
from mpl_toolkits.mplot3d import Axes3D


# Parameters for the plot.
iz_state = 70
x_skip = 5
y_skip = 50
z_skip = 5
cl = 10
r_min = 0.0
r_max = 0.1


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
#par_locs = par_locs[np.where(par_locs>0)]
#par_locs = par_locs[np.where(par_locs<cl)]

# Also mirror the par_locs to cover both sides (i.e. from -L to L instead of 0 to L).
# Need to add a zero in the as the middle point, hence two appends.
par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

# Get the DDLIM3 data to be used for the alpha variable. Normalize it.
print('Indexing DDLIM3...')
ddlim3 = netcdf.variables['DDLIM3'][::x_skip, iz_state, ::y_skip, ::z_skip].data
ddlim3  = ddlim3 / ddlim3.max()

fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')

# Use skip values to avoid overloading the plot and freezing the terminal.
x = pol_locs[::x_skip]
y = par_locs[::y_skip]
z = rad_locs[::z_skip]
print('Creating mesh...')
X, Y, Z = np.meshgrid(x, y, z)
points = np.array((X.flatten(), Y.flatten(), Z.flatten(), ddlim3.flatten()))
print("Creating plot...")
for p in range(0, len(points[0])):
    ax.scatter(points[0][p], points[1][p], points[2][p], color='r', alpha=points[3][p])
#ax.scatter(X, Y, Z)
#ax.set_ylim([-cl, cl])
#fig.tight_layout()
fig.show()

sum_range = np.where(np.logical_and(rad_locs > r_min, rad_locs < r_max))
