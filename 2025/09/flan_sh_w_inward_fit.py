import flan_plots
import matplotlib.pyplot as plt
import numpy as np


# Load case
path = "/home/zamp/flandir/sh_v6/sh_v6_w_inward.nc"
fp = flan_plots.FlanPlots(path)

# Average radial density
tmin_idx = -30
tmax_idx = -1
zidx = 8
t = fp.nc["time"][:].data
x = fp.nc["x"][:].data
nz = fp.nc["nz"][tmin_idx:tmax_idx, :, :, zidx].data.mean(axis=0).mean(axis=1)

# Fit range
fit_xmin = 2.31
fit_xmax = 2.345
fit_xmin_idx = np.where(x == x[x >= fit_xmin][0])[0][0]
fit_xmax_idx = np.where(x == x[x <= fit_xmax][-1])[0][0]
x_subset = x[fit_xmin_idx:fit_xmax_idx]
nz_subset = nz[fit_xmin_idx:fit_xmax_idx]

# Linear fit to log for decay rate
ln_nz_subset = np.log(nz_subset)
m, b = np.polyfit(x_subset, ln_nz_subset, 1)

# Data for plot
x_fit = np.linspace(x[fit_xmin_idx], x[fit_xmax_idx], 100)
nz_fit = np.exp(b) * np.exp(x_fit * m)

print("Decay length = {:.2f} cm".format(1/m * 100))

# Publication-ready plot
rsep = 2.259
lw = 3
fontsize = 16
fig, ax1 = plt.subplots(figsize=(6, 5))
ax1.plot(x-rsep, nz, color="tab:red", lw=lw)
ax1.plot(x_fit-rsep, nz_fit, color="k", linestyle="--", lw=lw)
ax1.set_yscale("log")
ax1.set_ylabel("W Density (arb.)", fontsize=fontsize)
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.tick_params(which="both", labelsize=fontsize-2)
ax1.grid(alpha=0.5)
fig.tight_layout()
fig.show()

