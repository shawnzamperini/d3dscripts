import numpy as np
import matplotlib.pyplot as plt
import netCDF4


# Load netcdf file.
#ncpath = "/Users/zamperini/Documents/lim_runs/bumper-001.nc"
ncpath = "/Users/zamperini/Documents/lim_runs/colprobe-a8.nc"
nc = netCDF4.Dataset(ncpath)

#def plot_rad_par_cont(plot_data, rmin=-0.01, rmax=0.00):
"""
Radial vs. parallel contour plot of the selected data.

plot_data (str): One of "ne", "te", "nz".
"""

plot_data = "nz"
rmin = 0.00
rmax = 0.02

# Load in the coordinates.
xs = nc.variables['XS'][:].data
xwids = nc.variables['XWIDS'][:].data
rad_locs = xs - xwids / 2.0
ps = nc.variables['PS'][:].data
pwids = nc.variables['PWIDS'][:].data
pol_locs = ps - pwids / 2.0
ys = nc.variables['YS'][:].data
ywids = nc.variables['YWIDS'][:].data
par_locs = ys - ywids / 2.0

# Mirror the parallel coordinate bc 3DLIM sucks.
par_locs = np.append(np.append(-par_locs[::-1], 0), par_locs)

# Impurity density.
if plot_data.lower() == "nz":
    ddlim3 = nc.variables['DDLIM3'][:].data.sum(axis=1)

    # Sum over the radial range to create a 2D plot.
    sum_range = np.where(np.logical_and(rad_locs > rmin, rad_locs < rmax))[0]
    summed_ddlim3 = ddlim3[:, :, sum_range].sum(axis=2)
    Z = summed_ddlim3

if plot_data.lower() == "ne":
    ne = nc.variables["CRNBS"][:]
    tmp_Z = np.tile(ne, 63).reshape((len(pol_locs), len(par_locs), len(rad_locs)))


# Plotting commands.
X, Y = np.meshgrid(par_locs, pol_locs)

fig, ax = plt.subplots()
#ax.contourf(X, Y, Z)
ax.pcolormesh(X, Y, Z)
ax.set_ylim([-0.1, 0.1])
fig.tight_layout()
fig.show()
