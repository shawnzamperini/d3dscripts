import netCDF4
import numpy as np
import pretty_plots as pp


test = '44'
filename = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8_test'+test+'.nc'
net = netCDF4.Dataset(filename)
print('Test: ' + test)

# 2D array of deposition data.
dep_arr = np.array(net.variables['NERODS3'][0] * -1)

# Location of each P bin, and its width. Currently they all have the same width,
# but it may end up such that there are custom widths so we leave it like this.
ps     = np.array(net.variables['PS'][:].data)
pwids  = np.array(net.variables['PWIDS'][:].data)

# Array of poloidal locations (i.e. the center of each P bin).
pol_locs = ps - pwids/2.0

# Distance cell centers along surface (i.e. the radial locations).
rad_locs = np.array(net.variables['ODOUTS'][:].data)

# Drop last row since it's garbage.
dep_arr = dep_arr[:-1, :]
pol_locs = pol_locs[:-1]

# Get data across probe closest to R = 10mm.
pol_slice = dep_arr[:, np.where(rad_locs > 0)[0][0]]

fig = pp.pplot(pol_locs, pol_slice, xrange=[-0.02,0.02])
