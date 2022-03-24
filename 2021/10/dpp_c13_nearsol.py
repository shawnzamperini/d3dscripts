# Plot the near-SOL C13 distribution for each flow case.
# Simple script to plot far-SOL distributions of C13 just to get a sense of
# what the C13 distribution feeding into the SOL could look like.
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge")
import oedge_plots
from scipy.signal import savgol_filter


include_lim_dist = False
lim_l = 13.79 + 32.39
source_start = 13.79 - 2.00
source_end = 13.79 - 0.01
#bkg_subs = [9.0e14, 6.5e14, 6.5e14, 6.5e14, 6.5e14]
bkg_subs = [0, 0, 0, 0, 0]

farsol = 13  # Near the edge of the grid.
smooth = False
ncpaths = [
"/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-007.nc",  # No flows
"/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-009_2.nc",  # M = 0.1
"/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-010_2.nc",  # M = 0.2
"/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-008_2.nc",  # M = 0.3
"/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-011_2.nc"   # M = 0.4
]
flows = [0.0, 0.1, 0.2, 0.3, 0.4]
ops = [oedge_plots.OedgePlots(nc) for nc in ncpaths]
ss = []
c_fars = []

for i in range(0, len(ops)):
    s, c_far = ops[i].along_ring(farsol, "DDLIMS", charge="all", plot_it=False)
    ss.append(s)
    c_fars.append(c_far)

# Scale the 3DLIM data onto this flux tube, background subtract the DIVIMP
# results, normalize.
max_c = 0
if include_lim_dist:
    lim_x = np.linspace(0, lim_l, 1000)
    lim_y = np.zeros(lim_x.shape)
    idx = np.where(np.logical_and(lim_x >= source_start, lim_x <= source_end))[0]
    lim_y[idx] = 1.0  # Rectangular injection source.
    lim_area = np.trapz(lim_y, lim_x)
    #lim_y = lim_y / lim_area  # Normalize.
    # Scale to DIVIMP ring size.
    lim_x = lim_x * max(s) / lim_x.max()

    for i in range(0, len(ops)):
        ss[i] = np.array(ss[i])
        keep = ss[i] <= 40
        c_fars[i] = np.array(c_fars[i])[keep] - bkg_subs[i]
        ss[i] = ss[i][keep]  # Remove data beyond where it gets near IT.

        if smooth:
            c_fars[i] = savgol_filter(c_fars[i], 11, 2)

        # Normalize to max between [5, 40].
        tmp = c_fars[i][ss[i] >= 5].max()
        if tmp > max_c:
            max_c = tmp
    for i in range(0, len(ops)):
        c_fars[i] = c_fars[i] / max_c

lw = 3
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, len(ops)))

fig, ax = plt.subplots()

for i in range(0, len(ss)):
    if smooth:
        if include_lim_dist:
            # Already smoothed above so we could normalize correctly.
            y = c_fars[i]
        else:
            y = savgol_filter(c_fars[i], 11, 2)
    else:
        y = c_fars[i]
    ax.plot(ss[i], y,
      label="M = {:.1f}".format(flows[i]), c=colors[i], lw=lw)

ax.legend(fontsize=14)
ax.set_xlabel("Distance from outer target (m)", fontsize=14)

# OMP location at about S = 0.5 / 2.0 of the total ring length.
s_omp = ss[i].max() * 0.5 / 2.0
s_crown = ss[i].max() / 2.0
ax.axvline(s_omp, color="k", linestyle="--")
ax.axvline(s_crown, color="k", linestyle="--")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.set_xlim([0, 86])
if include_lim_dist:
    ax.set_ylim([0, 1.1])
    #ax.plot(lim_x, lim_y)
    ax.fill_between(lim_x, lim_y, np.zeros(lim_y.shape), color="grey")
    ylabel = r"Normalized $\mathdefault{^{13}C\ Density\ }$"
else:
    ax.set_ylim([0.5e15, 6.0e15])
    #ax.set_yticks(np.arange(0.5e15, 3.0e15, 0.5e15))
    ylabel = r"$\mathdefault{^{13}C\ Density\ m^{-3}}$"
ax.set_ylabel(ylabel, fontsize=14)
#ax.set_yticklabels(np.arange(0.5e15, 3.0e15, 0.5e15))
#ax.set_yscale("log")
fig.tight_layout()
fig.show()
