# A script to show how to parse 3DLIM data to pull out the simulated deposition
# patterns.
import numpy as np
import matplotlib.pyplot as plt
import netCDF4


# Input parameters and options.
#ncpath = "/Users/zamperini/Documents/3DLIM Tutorial/lim_tutorial_5.nc"
#ncpath = "lim_tutorial_1.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_36.nc"
log_yscale = False

# Load netcdf file. Can investigate what is in it in the python interpretor
# with nc.variables.keys()
nc = netCDF4.Dataset(ncpath)

# Location of each P bin, and its width. Note syntax in pulling out data.
ps     = nc.variables['PS'][:].data
pwids  = nc.variables['PWIDS'][:].data

# Array of poloidal locations (i.e. the center of each P bin).
pol_locs = ps - pwids / 2.0

# Distance cell centers along probe surface (i.e. the radial locations).
rad_locs = nc.variables['ODOUTS'][:].data * 100  # m to cm

# Get the centerline index (or closest to it).
cline = np.abs(pol_locs).min()

# This is the deposition array of the 2D collector probe faces.
dep_arr = nc.variables['NERODS3'][0] * -1

# Index the deposition array at the centerline for plotting.
side1_x = rad_locs[np.where(rad_locs > 0.0)[0]]
side1_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
side2_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
side2_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

# Plotting commands.
fig, ax = plt.subplots()
ax.plot(side1_x, side1_y, lw=2, color="tab:red", label="Side 2")
ax.plot(side2_x, side2_y, lw=2, color="tab:purple", label="Side 1")
ax.legend(fontsize=14)
ax.set_xlabel("Distance along probe (cm)", fontsize=14)
ax.set_ylabel("Deposition (arbitrary units)", fontsize=14)
ax.grid()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
#ax.set_ylim([0.001, 1])
if log_yscale:
    ax.set_yscale("log")
fig.tight_layout()
fig.show()
