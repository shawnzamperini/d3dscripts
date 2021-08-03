# This script will generate a grid of plots showing key parameters for the
# two similar shots.
from gadata import gadata
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d


tags = ["BT0", "DENSV2", "PINJ", "PRAD_CORE", "FS04", "IP", "RVSOUT"]
signals247 = {}
signals277 = {}
for tag in tags:
    ga247 = gadata(tag, 167247)
    ga277 = gadata(tag, 167277)
    signals247[tag] = [ga247.xdata, ga247.zdata]
    signals277[tag] = [ga277.xdata, ga277.zdata]

# Fix weird shit with IP.
signals247["IP"] = [signals247["IP"][0].value, signals247["IP"][1].value]
signals277["IP"] = [signals277["IP"][0].value, signals277["IP"][1].value]

# Smooth PINJ.
signals247["PINJ"][1] = savgol_filter(signals247["PINJ"][1], 51, 2)
signals277["PINJ"][1] = savgol_filter(signals277["PINJ"][1], 51, 2)

# Calculate PSOL.
f_pinj247 = interp1d(signals247["PINJ"][0], signals247["PINJ"][1], bounds_error=False, fill_value=0)
f_pinj277 = interp1d(signals277["PINJ"][0], signals277["PINJ"][1], bounds_error=False, fill_value=0)
psol247 = f_pinj247(signals247["PRAD_CORE"][0]) * 1000 - signals247["PRAD_CORE"][1]
psol277 = f_pinj277(signals277["PRAD_CORE"][0]) * 1000 - signals277["PRAD_CORE"][1]
signals247["PSOL"] = [signals247["PRAD_CORE"][0], psol247]
signals277["PSOL"] = [signals277["PRAD_CORE"][0], psol277]

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 3

# Make a grid of plots with the signals. The bottom will be zoomed in on FS04.
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(5, 7))
fig, (ax2, ax3, ax4, ax5) = plt.subplots(4, 1, figsize=(5, 7))
#fig.subplots_adjust(hspace=0.25)

ax1.plot(signals277["BT0"][0], signals277["BT0"][1], color=colors[2], lw=lw, label="167277")
ax1.plot(signals247["BT0"][0], signals247["BT0"][1], color=colors[3], lw=lw, label="167247")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_ylabel(r"$\mathdefault{B_T}$ (T)", fontsize=fontsize)
ax1.set_xlim([0, 6000])
ax1.grid()
ax1.set_ylim([-3, 3])
ax1.legend(fontsize=12, loc="center right")
#ax1.tick_params(axis="x", which="both", labelbottom=False)
ax1.text(0.025, 0.8, "a)", fontsize=fontsize, transform=ax1.transAxes)

ax2.plot(signals277["DENSV2"][0], signals277["DENSV2"][1], color=colors[2], lw=lw)
ax2.plot(signals247["DENSV2"][0], signals247["DENSV2"][1], color=colors[3], lw=lw)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.set_ylabel(r"$\bar{\mathdefault{n}}_{\mathdefault{e}}\ (\mathdefault{m}^{-\mathdefault{3}})$", fontsize=fontsize)
ax2.set_xlim([0, 6000])
ax2.grid()
#ax2.tick_params(axis="x", which="both", labelbottom=False)
ax2.text(0.025, 0.8, "a)", fontsize=fontsize, transform=ax2.transAxes)

#ax3.plot(signals247["PINJ"][0], signals247["PINJ"][1]/1e3, color=colors[2], lw=lw)
#ax3.plot(signals277["PINJ"][0], signals277["PINJ"][1]/1e3, color=colors[3], lw=lw)
#ax3.spines["top"].set_visible(False)
#ax3.spines["right"].set_visible(False)
#ax3.set_ylabel(r"$\mathdefault{P_{inj}}$ (MW)", fontsize=fontsize)

ax3.plot(signals277["PSOL"][0], signals277["PSOL"][1]/1e6, color=colors[2], lw=lw)
ax3.plot(signals247["PSOL"][0], signals247["PSOL"][1]/1e6, color=colors[3], lw=lw)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.set_ylabel(r"$\mathdefault{P_{SOL}}$ (MW)", fontsize=fontsize)
ax3.set_xlim([0, 6000])
ax3.grid()
ax3.set_yticks([0, 3, 6])
#ax3.tick_params(axis="x", which="both", labelbottom=False)
ax3.text(0.025, 0.8, "b)", fontsize=fontsize, transform=ax3.transAxes)

ax4.plot(signals277["RVSOUT"][0], signals277["RVSOUT"][1], color=colors[2], lw=lw)
ax4.plot(signals247["RVSOUT"][0], signals247["RVSOUT"][1], color=colors[3], lw=lw)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)
ax4.set_ylabel("Strike\nPoint (m)", fontsize=fontsize)
ax4.set_xlim([0, 6000])
ax4.set_ylim(1.25, 1.45)
tile_x = [-1, 8000]
tile_l = np.full(len(tile_x), 1.32)
tile_r = np.full(len(tile_x), 1.37)
ax4.fill_between(tile_x, tile_l, tile_r, alpha=0.3, color=colors[1], edgecolor="none")
ax4.grid()
ax4.text(3500, 1.4, "W Ring", fontsize=fontsize, bbox=dict(facecolor="white", edgecolor="none"), color=colors[1])
#ax4.set_xlabel("Time (ms)", fontsize=fontsize)
ax4.text(0.025, 0.8, "c)", fontsize=fontsize, transform=ax4.transAxes)

ax5.plot(signals277["FS04"][0], signals277["FS04"][1], color=colors[2], lw=lw)
ax5.plot(signals247["FS04"][0], signals247["FS04"][1], color=colors[3], lw=lw)
ax5.spines["top"].set_visible(False)
ax5.spines["right"].set_visible(False)
ax5.set_ylabel(r"$\mathdefault{D}_\alpha$" + "\n" + r"(ph/$\mathdefault{cm^2}$/sr/s)", fontsize=fontsize)
ax5.set_xlim([2700, 2900])
ax5.grid()
ax5.set_ylim([0, 1.2e17])
ax5.set_xlabel("Time (ms)", fontsize=fontsize)
ax5.text(0.025, 0.8, "d)", fontsize=fontsize, transform=ax5.transAxes)

fig.tight_layout()
fig.show()
