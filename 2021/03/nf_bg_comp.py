# This script will make a plot of the along ring comparisons for ne and Te to
# Thomson scattering data. Also do a comparison of the radial profiles at the
# location the TS points mapped to the OMP.
# 2D grid where the left two are radial comparisons and the right two are along
# ring for ring 40.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Get background data.
unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
unf = oedge_plots.OedgePlots(unf_path)

# Load Thomson comparison file.
ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/ts167247_final_v2.xlsx"
ts_df = pd.read_excel(ts_path)
core_ts = ts_df[ts_df["System"] == "core"]

# Add the ne and Te values in the cell that the measurement falls under.
tes = []; nes = []; omps = []
for i in range(0, len(core_ts)):

    # Subtract 1 for indexing.
    ring = core_ts.iloc[i]["Ring"] - 1
    knot = core_ts.iloc[i]["Cell"] - 1

    # Grab the Te and ne values, as well as the midplane distance.
    tes.append(unf.nc.variables["KTEBS"][:][ring][knot])
    nes.append(unf.nc.variables["KNBS"][:][ring][knot])
    omps.append(unf.nc.variables["MIDIST"][:][1][ring])

# Put back into the dataframe.
core_ts["OSM Te"] = tes
core_ts["OSM ne"] = nes
core_ts["RMRSOMP"] = omps

# Can maybe get away with things by just sorting by RMRSOMP.
core_ts = core_ts.sort_values("RMRSOMP")

# Remove all core data.
core_ts = core_ts[core_ts["RMRSOMP"] > 0]

# Grab along ring data for ring 40 for both systems.
ring = 41
div_ts = ts_df[ts_df["System"] == "divertor"]
ring_core = core_ts[core_ts["Ring"] == ring - 1]
ring_div  = div_ts[div_ts["Ring"] == ring - 1]
osm_s = unf.nc.variables["KSS"][:][ring - 1]
osm_te = unf.nc.variables["KTEBS"][:][ring - 1]
osm_ne = unf.nc.variables["KNBS"][:][ring - 1]

# Drop zeros.
mask = osm_s != 0
osm_s = osm_s[mask]
osm_te = osm_te[mask]
osm_ne = osm_ne[mask]

# Average and std_dev of TS values.
core_avg_s = ring_core["S (m)"].mean()
core_avg_te = ring_core["Te (eV)"].mean()
core_std_te = ring_core["Te (eV)"].std()
core_avg_ne = ring_core["ne (m-3)"].mean()
core_std_ne = ring_core["ne (m-3)"].std()
div_avg_s = ring_div["S (m)"].mean()
div_avg_te = ring_div["Te (eV)"].mean()
div_std_te = ring_div["Te (eV)"].std()
div_avg_ne = ring_div["ne (m-3)"].mean()
div_std_ne = ring_div["ne (m-3)"].std()

# Plot it up.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 16
lw = 5
s = 100

# Plot the radial comparisons.
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,5))

ax1.scatter(core_ts["RMRSOMP"] * 100, core_ts["Te (eV)"], s=8, alpha=0.4,
  color=colors[2], marker="^")
ax1.plot(core_ts["RMRSOMP"] * 100, core_ts["OSM Te"], lw=lw, color=colors[1])
ax1.set_xlim([0, 7.5])
ax1.set_ylim([0, 80])
ax1.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax1.grid()
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(axis="x", which="both", labelbottom=False)
ax1.tick_params(axis="y", which="both", labelsize=13)
ax1.set_yticks(np.arange(0, 80, 20))
legend_elements = [Line2D([0], [0], color=colors[1], lw=lw, label="OSM"),
                   Line2D([0], [0], color=colors[2], marker="^", markersize=8, lw=0, label="Thomson")]
ax1.legend(handles=legend_elements, fontsize=13)

ax3.scatter(core_ts["RMRSOMP"] * 100, core_ts["ne (m-3)"], s=8, alpha=0.4,
  color=colors[2], marker="^")
ax3.plot(core_ts["RMRSOMP"] * 100, core_ts["OSM ne"], lw=lw, color=colors[1])
ax3.set_xlim([0, 7.5])
ax3.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (cm)", fontsize=fontsize)
ax3.set_ylabel(r"$\mathdefault{n_e\ (m}^{-\mathdefault{3}})$", fontsize=fontsize)
ax3.grid()
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.tick_params(axis="both", which="both", labelsize=13)

ax2.plot(osm_s, osm_te, color=colors[1], lw=lw)
ax2.scatter(ring_core["S (m)"], ring_core["Te (eV)"], marker=".", s=s,
  color=colors[2], alpha=0.3, zorder=50)
ax2.errorbar(core_avg_s, core_avg_te, yerr=core_std_te, marker=".", ms=11,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[2], zorder=51)
ax2.scatter(ring_div["S (m)"], ring_div["Te (eV)"], marker="o", s=s,
  facecolors='none', edgecolors=colors[2], alpha=0.3, zorder=52)
ax2.errorbar(div_avg_s, div_avg_te, yerr=div_std_te, marker=".", ms=11,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[2], zorder=53)
ax2.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax2.grid()
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(axis="x", which="both", labelbottom=False)
ax2.tick_params(axis="y", which="both", labelsize=13)
legend_elements = [Line2D([0], [0], color=colors[1], lw=lw, label="OSM"),
                   Line2D([0], [0], color=colors[2], marker=".", markersize=8, lw=0, label="Core"),
                   Line2D([0], [0], color=colors[2], marker="o", markersize=8, mfc="none", lw=0, label="Divertor")]
ax2.legend(handles=legend_elements, fontsize=13, loc="upper left")

ax4.scatter(ring_core["S (m)"], ring_core["ne (m-3)"], marker=".", s=s,
  color=colors[2], alpha=0.3, zorder=50)
ax4.errorbar(core_avg_s, core_avg_ne, yerr=core_std_ne, marker=".", ms=11,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[2], zorder=51)
ax4.scatter(ring_div["S (m)"], ring_div["ne (m-3)"], marker="o", s=s,
  facecolors='none', edgecolors=colors[2], alpha=0.3, zorder=52)
ax4.errorbar(div_avg_s, div_avg_ne, yerr=div_std_ne, marker=".", ms=11,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[2], zorder=53)
ax4.plot(osm_s, osm_ne, color=colors[1], lw=lw)
ax4.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax4.set_ylabel(r"$\mathdefault{n_e\ (m}^{-\mathdefault{3}})$", fontsize=fontsize)
ax4.grid()
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)
ax4.tick_params(axis="both", which="both", labelsize=13)

fig.tight_layout()
fig.show()
