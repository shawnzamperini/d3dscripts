# This script is just a test bed to see how to go from a 3DLIM simulation to a
# slice in the R, Z plane. The slice can only be taken at the origin I think.
# I can mostly justfiy that, but I can't really wrap my ahead around moving
# away from the origin since once you do that (i.e. Y != 0), you'd need to use
# updated connection lengths at whatever toroidal angle you'd then be at.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from collections import Counter
from matplotlib.colors import LogNorm


ncpath = "/Users/zamperini/Documents/lim_runs/actual-184527-004.nc"
r_tip = 2.29   # The R value of the origin to convert to easier coordinates.
p_to_z = [0.97636342, -0.188]  # Z = [0] * pbin + [1], i.e. linear fit to go
                               # from P to Z as figured out in
                               # bounds_file_from_mafot.py.

# Load the NetCDF and pull out various arrays.
print("Loading netCDF file and extracting variables...")
nc = netCDF4.Dataset(ncpath)
ps = nc.variables["PS"][:].data
xs = nc.variables["XS"][:].data
ys = nc.variables["YS"][:].data
ywids = nc.variables["YWIDS"][:].data
xwids = nc.variables["XWIDS"][:].data
pwids = nc.variables["PWIDS"][:].data
ddlim3 = nc.variables["DDLIM3"][:].data
vp = nc.variables["velplasma_4d_1"][:].data

# Warning if the width of the sum range is smaller than the smallest radial
# bin width...
min_xwid = xwids[np.nonzero(xwids)].min()

# Sum across all charge states.
print("Summing across charge states...")
ddlim3 = ddlim3.sum(axis=1)

# Calculate centers of bins. Special treatment for the Y coordinate since it
# needs to be mirrored and a zero added.
rad_locs = r_tip - (xs - xwids / 2)
pol_locs = p_to_z[0] * (ps - pwids / 2) + p_to_z[1]  # This is actually Z now.
tmp = ys - ywids / 2
par_locs = np.append(np.append(-tmp[::-1], 0), tmp)

y_keep_start = np.nonzero(par_locs)[0].min()
y_keep_end   = np.nonzero(par_locs)[0].max() + 1
x_keep_start = np.nonzero(rad_locs)[0].min()
x_keep_end   = np.nonzero(rad_locs)[0].max() + 1
p_keep_start = np.nonzero(pol_locs)[0].min()
p_keep_end   = np.nonzero(pol_locs)[0].max() + 1
rad_locs = rad_locs[x_keep_start:x_keep_end]
pol_locs = pol_locs[p_keep_start:p_keep_end]
par_locs = par_locs[y_keep_start:y_keep_end]
ddlim3 = ddlim3[p_keep_start:p_keep_end, y_keep_start:y_keep_end,
  x_keep_start:x_keep_end]

# Need to trim down rad_locs to get rid of the repeating values.
c = Counter(); idx = 0
for i in rad_locs:
    c[i] += 1
    if c[i] > 1:
        break
    idx += 1
new_rad_locs = rad_locs[:idx-1]

# Here is the new stuff. For density at the origin we will just do the average
# of the origin bin and the two surrounding it.
mid = int(len(par_locs) / 2)
ddlim3_mid = ddlim3[:, mid-1:mid+2, :idx-1].sum(axis=1)  # Shape is (pol_loc, rad_locs).
#ddlim3_mid = ddlim3[:, mid+1, :idx-1]

fig, ax = plt.subplots()
ax.pcolormesh(new_rad_locs, pol_locs, ddlim3_mid, shading="auto",
  norm=LogNorm(), cmap="inferno")
ax.set_xlabel("R (m)")
ax.set_ylabel("Z (m)")
#ax.set_aspect("equal")
fig.tight_layout()
fig.show()
