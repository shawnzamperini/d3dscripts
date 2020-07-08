import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# Some constants.
rsep = 2.20   # R of separatrix at Z = -0.188.

# Path to a file with the background data. This has 011c data.
path = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/d3d-167247-inj-021c.nc'
op = oedge_plots.OedgePlots(path)

# Get the radial Te, ne data at the probe location.
r, te = op.fake_probe(2.19, 2.3, -0.188, -0.188, 'Te', plot='R', show_plot=False)
r, ne = op.fake_probe(2.19, 2.3, -0.188, -0.188, 'ne', plot='R', show_plot=False)

# Make sure they're numpy arrays.
r  = np.array(r)
te = np.array(te)
ne = np.array(ne)

# Indices of that data beyond the separatrix and is nonzero.
te_idx = np.where(np.logical_and(r>rsep, te>1))

# Specific ne idx's to drop.
ne_idx = np.delete(te_idx, [0,3,4,5,6,7,8,9])

r_te = r[te_idx]
r_ne = r[ne_idx]
te   = te[te_idx]
ne   = ne[ne_idx]

# Use R-Rsep.
rminrsep_te = r_te - rsep
rminrsep_ne = r_ne - rsep

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

# Exponential fits to the data.
popt_te, pcov_te = curve_fit(exp_fit, rminrsep_te, te)
popt_ne, pcov_ne = curve_fit(exp_fit, rminrsep_ne, ne/1e19)
rminrsep_fit = np.linspace(0, 0.3, 100)
te_fit = exp_fit(rminrsep_fit, *popt_te)
ne_fit = exp_fit(rminrsep_fit, *popt_ne) * 1e19

# Output to be copy/pasted into 3DLIM input file.
zero_shift = 0.08
r_lim = zero_shift - rminrsep_fit

print("Te data")
for i in range(0, len(r_lim)):
    print("{:.3f}  {:.2f}".format(r_lim[i], te_fit[i]))
print("\nne data")
for i in range(0, len(r_lim)):
    print("{:.3f}  {:.2e}".format(r_lim[i], ne_fit[i]))
# Plot it.
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(15,5))
ax1.plot(rminrsep_te, te, '.')
ax2.plot(rminrsep_ne, ne, '.')
ax1.plot(rminrsep_fit, te_fit, '-')
ax2.plot(rminrsep_fit, ne_fit, '-')
ax1.set_xlim([0, 0.15])
ax1.set_xlabel('R-Rsep (m)')
ax2.set_xlabel('R-Rsep (m)')
ax1.set_ylabel('Te (eV)')
ax2.set_ylabel('ne (m-3)')
fig.tight_layout()
fig.show()
