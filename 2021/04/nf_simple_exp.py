import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
import oedge_plots

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)


unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-035a.nc"
unf = oedge_plots.OedgePlots(unf_path)
unf_ring = 70
unf_along = list(unf.along_ring(unf_ring, "DDLIMS", charge="all", plot_it=False))
unf_along_raw = unf_along.copy()
unf_along[1] = savgol_filter(unf_along[1], 19, 2)
unf_along[1][unf_along[0] > 42] = unf_along_raw[1][unf_along[0] > 42]
unf_y = unf_along[1][np.logical_and(unf_along[0] > 35, unf_along[0] < 42)]

x = np.linspace(0, 9)
y = np.exp(-x)
x_unf = np.linspace(0, 7, len(unf_y))

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))

fig, ax = plt.subplots(figsize=(6, 2))
ax.plot(x, y/y.max(), lw=3, color=colors[2])
#ax.plot(x_unf, (unf_y-unf_y.min())/(unf_y.max()-unf_y.min()))
ax.plot(x_unf, unf_y/unf_y.max(), color=colors[2], lw=3, linestyle="--")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.tick_params(axis="both", which="both", labelleft=False, labelbottom=False, left=False, bottom=False)
ax.set_yticks([0, 0.5, 1])
fig.tight_layout()
fig.show()
