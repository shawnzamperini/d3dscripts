# Purpose of this script is to plot the radial W density at various times
# to see when we can consider it steady-state.
import flan_plots
import matplotlib.pyplot as plt
import numpy as np


# Load the "normal" case, essentially the poster child plot. Variance reduction
# was on.
path = "/Users/zamperini/flandir/coll_on/coll_on.nc"
fp = flan_plots.FlanPlots(path)

# Lists to hold the resulting profiles.
nzs = []

frames = [300, 400, 500, 600]
for i in range(len(frames)):
    print(frames[i])
    fp_data = fp.plot_profiles(["vx", "nz"], plot_z=0.3125,
        normtype=["log","log"], vmin=[1e18, 1], vmax=[2e19, 1e3], 
        skip_video=True, show_plot=False, x_offset=-2.259, f=frames[i])
    
    # Get R-Rsep and nz. Only average non-zero cells.
    x = fp_data["x"]
    nz = fp_data["data"][1][frames[i]]
    nz[nz == 0] = np.nan
    nz = np.nanmean(nz, axis=1)

    nzs.append(nz)


fontsize = 14
fig, ax1 = plt.subplots(figsize=(5, 4))
for i in range(len(frames)):
    ax1.plot(x, nzs[i], label="{:.0f} us".format(frames[i] / 10))
ax1.legend()
ax1.set_yscale("log")
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=fontsize)
ax1.set_ylabel(r"W Density $\mathdefault{m^{-3}}$", fontsize=fontsize)
fig.tight_layout()
fig.show()
