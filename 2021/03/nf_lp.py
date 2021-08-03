# This script generates a plot of the Langmuir probe data and the fit to it
# that is input to OSM.
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Load the Excel file with the LP data and pull out the relvant data.
lp_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/lp_input_167247_v3.xlsx"
lp_df = pd.read_excel(lp_path, skiprows=1)

lp_psin  = lp_df["Psin.2"]
lp_te    = lp_df["Te (eV).3"]
lp_jsat  = lp_df["jsat (A/m2)"]
dts_psin = lp_df["Psin.3"]
dts_te   = lp_df["Te (eV).4"]
dts_jsat = lp_df["ne*cs/2"]
fit_psin = lp_df["Psin.1"]
fit_te   = lp_df["Te (eV).1"]
fit_jsat = lp_df["Jsat (A/m2).1"]

# Remove large dts data points.
mask     = dts_te < 10
dts_psin = dts_psin[mask]
dts_te   = dts_te[mask]
dts_jsat = dts_jsat[mask]

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 3
xlims = [0.98, 1.18]
s = 30

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.5))

ax1.plot(fit_psin, fit_te, color=colors[1], lw=lw, zorder=10)
ax1.scatter(lp_psin, lp_te, color=colors[2], marker="o", s=s, zorder=11, edgecolor="k")
ax1.scatter(dts_psin, dts_te, color=colors[3], marker="^", s=s, zorder=12, edgecolor="k")
ax1.set_xlim(xlims)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.grid()
ax1.set_xlabel("Psin", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax1.set_ylim([0, None])

ax2.plot(fit_psin, fit_jsat, color=colors[1], lw=lw, label="Input", zorder=10)
ax2.scatter(lp_psin, lp_jsat, color=colors[2], label="LP", marker="o", s=s, zorder=11, edgecolor="k")
ax2.scatter(dts_psin, dts_jsat, color=colors[3], label="DTS", marker="^", s=s, zorder=12, edgecolor="k")
ax2.set_xlim(xlims)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.grid()
ax2.legend(fontsize=12)
ax2.set_xlabel("Psin", fontsize=fontsize)
ax2.set_ylabel(r"$\mathdefault{j_{sat}\ (A/m^2)}$", fontsize=fontsize)
ax2.ticklabel_format(axis="y", style="scientific", scilimits=(4,4))
ax2.set_ylim([0, None])
ax2.set_yticks(np.arange(0, 210000, 30000))
ax2.text(1.17, 60000, "#167247", rotation=90, bbox=dict(facecolor="white", edgecolor="none"))

fig.tight_layout()
fig.show()
