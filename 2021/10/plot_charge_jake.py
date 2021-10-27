import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge")
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


ncpath = "/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-007.nc"
op = oedge_plots.OedgePlots(ncpath)

vals = []
for charge in range(0, 8):
    vals.append(op.along_ring(14, "DDLIMS", charge=charge, plot_it=False))

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 8))

fig, axs = plt.subplots(2, 4, figsize=(12,6))
axs = axs.flatten()
for i in range(0, len(vals)):
    if i == 0:
        label = "Primary Neutrals"
    elif i == 1:
        label = "Total Neutrals"
    else:
        label = "C{}+".format(i-1)
    axs[i].plot(vals[i][0], vals[i][1], label=label, color=colors[i], lw=2)
    axs[i].spines["top"].set_visible(False)
    axs[i].spines["right"].set_visible(False)
    axs[i].grid()
    axs[i].legend()
axs[4].set_xlabel("S (m)")
axs[0].set_ylabel("Density (m-3)")
#axs[0].legend()
#ax.set_yscale("log")
#ax.set_ylim(1e12, 1e16)
fig.tight_layout()
fig.show()
