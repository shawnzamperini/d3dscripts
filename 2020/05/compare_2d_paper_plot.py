import lim_plots as lim
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker, cm


# Some constants.
rad_cutoff = 0.05
probe_width = 0.015
vmax = 0.75

# Load the LimPlots objects for a simple and complex SOL case.
#simple_path  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-037.nc'
simple_path  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-046b.nc'
complex_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-036.nc'
simp = lim.LimPlots(simple_path)
comp = lim.LimPlots(complex_path)

# Grab the deposition arrays.
simp_dep = simp.get_dep_array()
comp_dep = comp.get_dep_array()

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
comp_dep = comp_dep[:, idx]

# Get only positive values of rad_locs for ITF...
idx = np.where(rad_locs > 0.0)[0]
X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
simp_Z_itf = simp_dep[:, idx]
comp_Z_itf = comp_dep[:, idx]

# ... negative for OTF.
idx = np.where(rad_locs < 0.0)[0]
X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
simp_Z_otf = simp_dep[:, idx][:, ::-1]
comp_Z_otf = comp_dep[:, idx][:, ::-1]

# Make the levels for the contour plot out of whichever side has the max deposition.
levels = np.linspace(0, vmax, 15)
#lev_exp = np.arange(0.1, 1)
#levels = np.power(10, lev_exp)

# Use the OTF data.
X = X_otf; Y = Y_otf

fig = plt.figure()

# Big plot for sharing the Y-axis. Turn off axis lines and ticks.
ax = fig.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax.set_ylabel('Poloidal (cm)\n', fontsize=16)

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.contourf(X*100, Y*100, simp_Z_otf/simp_Z_otf.max(), levels=levels, cmap='Reds')
ax2.contourf(X*100, Y*100, comp_Z_otf/comp_Z_otf.max(), levels=levels, cmap='Reds')
#ax1.contourf(X*100, Y*100, simp_Z_otf/simp_Z_otf.max(), locator=ticker.LogLocator(), cmap='Reds')
#ax2.contourf(X*100, Y*100, comp_Z_otf/comp_Z_otf.max(), locator=ticker.LogLocator(), cmap='Reds')
ax1.set_ylim([-1.5, 1.5])
ax2.set_ylim([-1.5, 1.5])

fig.tight_layout()
fig.show()
