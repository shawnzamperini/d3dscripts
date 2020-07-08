import lim_plots as lim
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker, cm


# Some constants.
rad_cutoff = 0.05
probe_width = 0.015
vmax = 0.75

# Load the LimPlots objects for a simple and complex SOL case.
simple_path  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-035.nc'
simp = lim.LimPlots(simple_path)

# Grab the deposition arrays.
simp_dep = simp.get_dep_array()

# Location of each P bin, and its width. Currently they all have the same width,
# but it may end up such that there are custom widths so we leave it like this.
ps     = np.array(simp.nc.variables['PS'][:].data)
pwids  = np.array(simp.nc.variables['PWIDS'][:].data)

# Array of poloidal locations (i.e. the center of each P bin).
pol_locs = ps - pwids/2.0

# Distance cell centers along surface (i.e. the radial locations).
rad_locs = np.array(simp.nc.variables['ODOUTS'][:].data)

# Remove data beyond rad_cutoff.
idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
rad_locs = rad_locs[idx]
simp_dep = simp_dep[:, idx]

# Get only positive values of rad_locs for ITF...
idx = np.where(rad_locs > 0.0)[0]
X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
simp_Z_itf = simp_dep[:, idx]

# ... negative for OTF.
idx = np.where(rad_locs < 0.0)[0]
X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
simp_Z_otf = simp_dep[:, idx][:, ::-1]

# Make the levels for the contour plot out of whichever side has the max deposition.
levels = np.linspace(0, vmax, 15)
#lev_exp = np.arange(0.1, 1)
#levels = np.power(10, lev_exp)

# Use the OTF data.
X = X_otf; Y = Y_otf

# Normalize.
max_val = max(simp_Z_itf.max(), simp_Z_otf.max())
simp_Z_itf = simp_Z_itf / max_val
simp_Z_otf = simp_Z_otf / max_val

otf_vmax = 0.25
simp_Z_otf = np.clip(simp_Z_otf, 0, otf_vmax)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
cf1 = ax1.contourf(X, Y, simp_Z_itf, cmap='Reds', vmin=0, vmax=1, levels=10)
cf2 = ax2.contourf(X, Y, simp_Z_otf, cmap='Purples', vmin=0, vmax=otf_vmax, levels=10)
fig.colorbar(cf1, ax=ax1)
fig.colorbar(cf2, ax=ax2)
ax1.set_ylim(-0.015, 0.015)
ax1.set_xlim(0, 0.05)
fig.tight_layout()
fig.show()
