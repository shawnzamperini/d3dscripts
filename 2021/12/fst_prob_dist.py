# Script for FST paper to show update in probability distribution.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

divimp_nc_path = "/Users/zamperini/Documents/d3d_work/167196/d3d-167196-mrc-shifted-nodrift.nc"
op = oedge_plots.OedgePlots(divimp_nc_path)

s, new = op.along_ring(91, "DDLIMS", charge="all", plot_it=False)

old = np.zeros(len(new))
old[np.logical_and(s>3, s<8)] = 1.0

# Normalize both.
new_sum = np.trapz(new, s)
old_sum = np.trapz(old, s)
new = new / new_sum
old = old / old_sum

# A flat distribution for comparison.
ones = np.ones(s.shape)
flat = ones / np.trapz(ones, s)

# Ticks.
xticks = np.arange(0, 25, 5)
yticks = np.arange(0, 0.25, 0.05)

fig, ax = plt.subplots(figsize=(5,4))
ax.plot(s, new, label="DIVIMP (New)", color="tab:green", lw=2, zorder=20)
ax.plot(s, old, label="Rectangular (Old)", color="tab:pink", lw=2, zorder=10)
ax.plot(s, flat, label="Flat", color="tab:orange", lw=2, zorder=5)
ax.legend(fontsize=12)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_ylim([0, None])
ax.set_xlabel("Distance from upper baffle (m)", fontsize=12)
ax.set_ylabel("Injection probability distribution", fontsize=12)
fig.tight_layout()
fig.show()
