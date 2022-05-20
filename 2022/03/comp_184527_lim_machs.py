# Script to compare the deposition profiles for a set of runs that scanned the
# parallel drift velocity from M = 0.2 (outer) to -0.1 (inner) for the 184527
# collector probe.
import LimPlots
import matplotlib.pyplot as plt
from importlib import reload
import pandas as pd
from scipy.signal import savgol_filter
import numpy as np


# This absfac is just the maximum of the input data for Z02.
absfac = 2.662e+15

# Load the data.
reload(LimPlots)
root = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/"
#case_nums = range(11, 18)
#machs = [0.20, 0.15, 0.10, 0.05, 0.00, -0.05, -0.10]
lps = []
for num in case_nums:
    path = root + "mcp-184527-0{}.nc".format(num)
    lps.append(LimPlots.LimPlots(path))

# Extract the deposition data. 184527: 1 = OTF, 2 = ITF
itf_xs = []; otf_xs = []
itf_ys = []; otf_ys = []
for i in range(0, len(lps)):
    cent = lps[i].centerline(showplot=False)
    itf_xs.append(cent["x2"] * 100)
    itf_ys.append(cent["y2"] * absfac)
    otf_xs.append(cent["x1"] * 100)
    otf_ys.append(cent["y1"] * absfac)

# Load in the LAMS data. Remove data before zero, smoothing, background subtraction.
xlpath = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"
df = pd.read_excel(xlpath, sheet_name="For Export")
lams_itf_x = df["ml04_loc"]
lams_itf_y = df["ml04_excess_c13"]
lams_otf_x = df["mr04_loc"]
lams_otf_y = df["mr04_excess_c13"]
mask1 = lams_itf_x > 0
mask2 = lams_otf_x > 0
lams_itf_x = lams_itf_x[mask1]
lams_itf_y = lams_itf_y[mask1]
lams_otf_x = lams_otf_x[mask2]
lams_otf_y = lams_otf_y[mask2]

# Convert LAMS counts to 1e17 atoms/cm2, and then to atoms/m2.
lams_itf_y = (lams_itf_y + 346.05) / 11942
lams_otf_y = (lams_otf_y + 346.05) / 11942
#lams_itf_y = lams_itf_y * 1e4 * 1e17
#lams_otf_y = lams_otf_y * 1e4 * 1e17
lams_itf_y = lams_itf_y * 1e17
lams_otf_y = lams_otf_y * 1e17

lams_itf_ys = savgol_filter(lams_itf_y, 51, 2)
lams_otf_ys = savgol_filter(lams_otf_y, 51, 2)
min_itfy = lams_itf_ys.min()
min_otfy = lams_otf_ys.min()
lams_itf_y = lams_itf_y - min_itfy
lams_otf_y = lams_otf_y - min_otfy
lams_itf_ys = lams_itf_ys - min_itfy
lams_otf_ys = lams_otf_ys - min_otfy



# Plotting.
cmap = plt.get_cmap("inferno")
colors = [cmap(j) for j in np.linspace(0, 1, len(lps))]
lim_ylim = [0, None]
lam_ylim = [0, None]
fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(10, 4))
ax2 = ax1.twinx()
ax4 = ax3.twinx()

ax1.plot(lams_itf_x, lams_itf_y, color="tab:red", alpha=0.3)
ax1.plot(lams_itf_x, lams_itf_ys, color="tab:red", label="ITF", lw=3)
ax3.plot(lams_otf_x, lams_otf_y, color="tab:purple", alpha=0.3)
ax3.plot(lams_otf_x, lams_otf_ys, color="tab:purple", label="OTF", lw=3)

for i in range(0, len(lps)):
    ax2.plot(itf_xs[i], itf_ys[i], color="k", lw=3)
    ax2.plot(itf_xs[i], itf_ys[i], color=colors[i], lw=2)
    ax4.plot(otf_xs[i], otf_ys[i], color="k", lw=3)
    ax4.plot(otf_xs[i], otf_ys[i], color=colors[i], lw=2, label="{:.2f}".format(machs[i]))

ax4.legend()
ax1.set_ylim(lam_ylim)
ax3.set_ylim(lam_ylim)
ax2.set_ylim(lim_ylim)
ax4.set_ylim(lim_ylim)
ax1.set_xlim([0, 10])
ax3.set_xlim([0, 10])
ax1.legend()
ax3.legend()
ax1.set_xlabel("Distance along probe (cm)")
ax3.set_xlabel("Distance along probe (cm)")
ax1.set_ylabel("LAMS (atoms/cm2)")
ax3.set_ylabel("LAMS (atoms/cm2)")
ax4.set_ylabel("3DLIM Counts")
fig.tight_layout()
fig.show()
