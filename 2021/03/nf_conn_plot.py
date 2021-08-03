import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Path to conneciton length data.
path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/Connection Lengths/167247/167247.xlsx"
path277 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/Connection Lengths/167277/167277.xlsx"
itf = pd.read_excel(path, sheet_name="MAFOT ITF", skiprows=2)
otf = pd.read_excel(path, sheet_name="MAFOT OTF", skiprows=2)
itf277 = pd.read_excel(path277, sheet_name="MAFOT ITF", skiprows=2)
otf277 = pd.read_excel(path277, sheet_name="MAFOT OTF", skiprows=2)

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 3

ytick_labels = [r"10$^{-\mathdefault{1}}$", r"10$^{\mathdefault{0}}$",
  r"10$^{\mathdefault{1}}$"]

fig, ax = plt.subplots(figsize=(5.5, 4))

# Color regions to align with the layout plot.
wall_patch = patches.Rectangle((0.25, 0), 0.75, 1, transform=ax.transAxes, color=colors[4], alpha=0.6)
ax.add_patch(wall_patch)
far_patch = patches.Rectangle((0, 0), 0.25, 1, transform=ax.transAxes, color=colors[3], alpha=0.6)
ax.add_patch(far_patch)

ax.plot(otf["R-Rsep OMP (cm)"], otf["Connection Length (km)"]*1000, label="OTF",
  lw=lw, color=colors[3])
ax.plot(otf277["R-Rsep OMP (cm)"], otf277["Connection Length (km)"]*1000,
  lw=lw, color=colors[3], linestyle="--")
ax.plot(itf["R-Rsep OMP (cm)"], itf["Connection Length (km)"]*1000, label="ITF",
  lw=lw, color=colors[2])
ax.plot(itf277["R-Rsep OMP (cm)"], itf277["Connection Length (km)"]*1000,
  lw=lw, color=colors[2], linestyle="--")
ax.set_yscale("log")
ax.set_xlim([5, 15])
ax.set_ylim([0.1, 40])
ax.legend(fontsize=13, loc=(0.02, 0.2))
ax.grid(which="both", alpha=0.6)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (cm)", fontsize=fontsize)
ax.set_ylabel("Connection Length (m)", fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_yticks([0.1, 1, 10])
ax.set_yticklabels(ytick_labels)
ax.text(6.3, 10.8, "Upper Divertor", rotation=0, fontsize=fontsize)
ax.text(7.9, 2.4, "Upper Baffle", rotation=-8, fontsize=fontsize)
ax.text(11, 0.7, "Outer Wall", rotation=-18, fontsize=fontsize)
ax.text(10.5, 8.5, "Outer Target", rotation=-5, fontsize=fontsize)
ax.text(11.5, 0.13, "#167247, 277", fontsize=fontsize)
ax.text(5.12, 0.13, "Far-SOL", fontsize=fontsize)
ax.text(8.3, 0.13, "Wall-SOL", fontsize=fontsize)

fig.tight_layout()
fig.show()
