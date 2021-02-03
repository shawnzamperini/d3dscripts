import matplotlib.pyplot as plt
import numpy as np
import oedge_plots
from scipy.optimize import curve_fit


# Load the divimp run.
ncpath1 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-025a.nc'
#ncpath2 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-025b.nc'
ncpath2 = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/utk-divimp/d3d-167247-inj-026d.nc'
op1 = oedge_plots.OedgePlots(ncpath1)
op2 = oedge_plots.OedgePlots(ncpath2)

# Get the along ring data of the W density.
s1, w1 = op1.along_ring(70, "DDLIMS", charge="all", plot_it=False)
s2, w2 = op2.along_ring(70, "DDLIMS", charge="all", plot_it=False)

# Restrict it to just the region below the upper baffle.
mask = s1 > 25
s1 = s1[mask]
s2 = s2[mask]
w1 = w1[mask][::-1]
w2 = w2[mask][::-1]

# Set s starting at 2.
s1 = s1 - s1.min()
s2 = s2 - s2.min()
mask = s1 > 2
s1 = s1[mask]
s2 = s2[mask]
w1 = w1[mask]
w2 = w2[mask]
s1 = s1 - s1.min()
s2 = s2 - s2.min()

# Exponential fit and all that.
def exp_fit(x, a, b, c):
    return a * np.exp(b * x) + c
#popt1, pcov = curve_fit(exp_fit, s1, w1)
#w1_fit = exp_fit(s1, *popt1)
#print("Lambda1 = {:.2f} m".format(1/popt1[1]))

# Make a plot to fit above the 3DLIM schematic on the slide.
fig, ax = plt.subplots(figsize=(7,3))
ax.plot(s1, w1, 'darkorange', lw=5)
#ax.plot(s1, w1_fit, 'k--', lw=3)
ax.plot(s2, w2, 'lightskyblue', lw=5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(np.arange(0, 18, 2))
ax.tick_params(which='both', labelsize=12)
ax.set_xlabel("Distance from outer target (m)", fontsize=14)
ax.set_ylabel("W Density\n(arbitrary units)", fontsize=14)
ax.set_yticks([])
fig.tight_layout()
fig.show()
