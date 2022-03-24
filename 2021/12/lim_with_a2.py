# Script to compare a 3DLIM run with A2.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4


# Input parameters and options.
ncpath = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_50a.nc"
log_yscale = False

a2_path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2_path)

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

a2_itf_x = a2["Distance from Tip D (cm)"].values
a2_otf_x = a2["Distance from Tip U (cm)"].values
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax3 = ax1.twiny()
ax4 = ax2.twiny()

ax4.plot(side1_x, side1_y, label="OTF (RBS)", zorder=50, color="tab:purple")
ax4.plot(side2_x, side2_y, label="ITF (RBS)", zorder=50, color="tab:red")
ax3.scatter(a2_itf_x, a2_itf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:red", label="ITF (3DLIM)")
ax3.scatter(a2_otf_x, a2_otf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple", label="OTF (3DLIM)")
ax4.legend()
ax4.set_ylim([0, 35])
ax1.set_xlabel("Distance along probe (cm)")
ax1.set_ylabel("W Areal Density (W/cm2)")
ax2.set_ylabel("3DLIM Deposition (arbitrary)")

ax3.tick_params(axis='x', which='both', bottom=False, top=False)
ax4.tick_params(axis='x', which='both', bottom=False, top=False)
ax3.set_xticklabels([])
ax4.set_xticklabels([])

ax1.set_zorder(1)
ax2.set_zorder(2)
ax3.set_zorder(3)
ax4.set_zorder(4)

fig.tight_layout()
fig.show()
