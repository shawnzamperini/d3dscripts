# This script generates a plot of the parallel profiles of the impurity density
# along the furthest far-SOL flux tube.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
import matplotlib as mpl
from scipy.signal import savgol_filter

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Get case representative of each direction.
#unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-035a.nc"
#fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031d.nc"
#fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-002c2.nc"

fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-008c.nc"

unf = oedge_plots.OedgePlots(unf_path)
fav = oedge_plots.OedgePlots(fav_path)

# Load along ring total impurity density data.
fav_ring = 42
#fav_ring = 70
unf_ring = 70
unf_along = list(unf.along_ring(unf_ring, "DDLIMS", charge="all", plot_it=False))
fav_along = list(fav.along_ring(fav_ring, "DDLIMS", charge="all", plot_it=False))

# Don't include all the junk at the inner target for favorable BT.
#mask = fav_along[0] > 10
#fav_along[0] = fav_along[0][mask]
#fav_along[1] = fav_along[1][mask]

# Get maximum impurity density value.
max_w = 0
for op in [unf_along, fav_along]:
    if op[1].data.max() > max_w:
        max_w = op[1].data.max()

# Normalize.
for op in [unf_along, fav_along]:
    #op[1] = op[1] / np.max(op[1])
    op[1] = op[1] / max_w

# Smooth away the artifact in the unf profile...
unf_along_raw = unf_along.copy()
unf_along[1] = savgol_filter(unf_along[1], 19, 2)
unf_along[1][unf_along[0] > 42] = unf_along_raw[1][unf_along[0] > 42]
fav_along_raw = fav_along.copy()
fav_along[1] = savgol_filter(fav_along[1], 19, 2)
fav_along[1][fav_along[0] > 42] = fav_along_raw[1][fav_along[0] > 42]

# Distance to separatrix at OMP.
unf_mid_dist = unf.nc.variables["MIDIST"][1][unf_ring]
unf_mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(unf_mid_dist*100)
fav_mid_dist = fav.nc.variables["MIDIST"][1][fav_ring]
fav_mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(fav_mid_dist*100)

# For the two rings 70 and 42 the distances are 6.86 and 6.62 cm, respectively.
# Just set an average of 6.7 cm.
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(6.7)

# Plot it up.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5

fig, ax = plt.subplots(figsize=(5.5, 4))
rect = patches.Rectangle((35, 0), 43.5-35, 1, color=colors[4], alpha=0.75)
ax.add_patch(rect)
ax.plot(unf_along[0], unf_along[1], label="Unfavorable", lw=lw, color=colors[2])
ax.plot(fav_along[0], fav_along[1], label="Favorable", lw=lw, color=colors[3])
ax.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax.set_ylabel("Tungsten density (normalized)", fontsize=fontsize)
ax.legend(fontsize=13, framealpha=1.0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.grid()
ax.set_xlim([15, None])
#ax.set_xticks(np.arange(25, 50, 5))
ax.set_ylim([0.0, 0.25])
#ax.axvline(33, color="k", linestyle="--", lw=lw)
ax.axvline(43.5, color="k", linestyle="--", lw=lw)
ax.text(45, 0.025, "W Injection Location", rotation=90, fontsize=fontsize, bbox=dict(color="white"))
#ax.text(30.5, 1.15, "Upper Baffle", rotation=0, fontsize=fontsize, bbox=dict(color="white"))
#ax.text(29.7, 1.05, "Above", rotation=0, fontsize=fontsize, bbox=dict(color="white"))
#ax.text(33.5, 1.05, "Below", rotation=0, fontsize=fontsize, bbox=dict(color="white"))
ax.text(16, 0.17, mid_str, fontsize=fontsize, bbox=dict(color="white"))
ax.text(33, 0.2, "Source near\nwall-SOL", rotation=0, fontsize=fontsize, bbox=dict(color=colors[4], ec="k"))

fig.tight_layout()
fig.show()
