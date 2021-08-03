# This script shows background comparisons but will also include the input
# LP data. The plot grid will be a 2x3.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


# Input to pick which shot.
shot = 167277

if shot == 167247:
    lp_jsat_mult = 1.25
    lp_te_mult = 1.0
    #op_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
    op_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a.nc"
    ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/ts167247_final_v2.xlsx"
    ts_ring = 41
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
    knot_mod = 0
    jsat_yticks = np.arange(0, 210000, 30000)
    ts_te_ylim = [None, 55]

elif shot == 167277:
    lp_jsat_mult = 1.25
    lp_te_mult = 1.40
    #op_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-bkg-003.nc"
    op_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
    ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/setup-files/ts167277_final_v3.xlsx"
    ts_ring = 31
    lp_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/setup-files/lp_167277.xlsx"
    lp_df = pd.read_excel(lp_path, sheet_name="Data", skiprows=1)
    lp_psin = lp_df["psin"]
    lp_te = lp_df["Te (eV)"]
    lp_jsat = lp_df["jsat (A/m2)"]
    dts_psin = lp_df["psin.1"]
    dts_te = lp_df["Te (eV).1"]
    dts_jsat = lp_df["ne*cs/2"]
    fit_psin = lp_df["psin.3"]
    fit_te = lp_df["Te"]
    fit_jsat = lp_df["jsat"]
    knot_mod = 0
    jsat_yticks = np.arange(0, 240000, 60000)
    ts_te_ylim = [None, 40]


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)
pd.options.mode.chained_assignment = None

# Get background data.
#unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
op = oedge_plots.OedgePlots(op_path)

# Load Thomson comparison file.
#ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/ts167247_final_v2.xlsx"
ts_df = pd.read_excel(ts_path)
core_ts = ts_df[ts_df["System"] == "core"]

# Add the ne and Te values in the cell that the measurement falls under.
tes = []; nes = []; omps = []
for i in range(0, len(core_ts)):

    # Subtract 1 for indexing.
    ring = core_ts.iloc[i]["Ring"] - 1
    knot = core_ts.iloc[i]["Cell"] - 1 + knot_mod

    # Grab the Te and ne values, as well as the midplane distance.
    tes.append(op.nc.variables["KTEBS"][:][ring][knot])
    nes.append(op.nc.variables["KNBS"][:][ring][knot])
    omps.append(op.nc.variables["MIDIST"][:][1][ring])

# Put back into the dataframe.
core_ts["OSM Te"] = tes
core_ts["OSM ne"] = nes
core_ts["RMRSOMP"] = omps

# Can maybe get away with things by just sorting by RMRSOMP.
core_ts = core_ts.sort_values("RMRSOMP")

# Remove all core data.
core_ts = core_ts[core_ts["RMRSOMP"] > 0]

# Grab along ring data for ring 40 for both systems.
#ring = 41
div_ts = ts_df[ts_df["System"] == "divertor"]
ring_core = core_ts[core_ts["Ring"] == ts_ring - 1]
ring_div  = div_ts[div_ts["Ring"] == ts_ring - 1]
osm_s = op.nc.variables["KSS"][:][ts_ring - 1]
osm_te = op.nc.variables["KTEBS"][:][ts_ring - 1]
osm_ne = op.nc.variables["KNBS"][:][ts_ring - 1]

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

# Distance from separatrix OMP.
mid_dist = op.nc.variables["MIDIST"][1][ts_ring]
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(mid_dist*100)
psin_val = op.nc.variables["PSIFL"][:][ts_ring][0]
psin_str = r"$\mathrm{\psi_N}$" + " = {:.2f}".format(psin_val)

# Load the Excel file with the LP data and pull out the relevant data.
#lp_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/lp_input_167247_v3.xlsx"
#lp_df = pd.read_excel(lp_path, skiprows=1)

#lp_psin  = lp_df["Psin.2"]
#lp_te    = lp_df["Te (eV).3"]
#lp_jsat  = lp_df["jsat (A/m2)"]
#dts_psin = lp_df["Psin.3"]
#dts_te   = lp_df["Te (eV).4"]
#dts_jsat = lp_df["ne*cs/2"]
#fit_psin = lp_df["Psin.1"]
#fit_te   = lp_df["Te (eV).1"]
#fit_jsat = lp_df["Jsat (A/m2).1"]

# Remove large dts data points.
mask     = dts_te < 10
dts_psin = dts_psin[mask]
dts_te   = dts_te[mask]
dts_jsat = dts_jsat[mask]

# Plot it up.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 12
lw = 5
s = 100
s_lp = 40
xlims = [0.98, 1.18]

# Plot the radial comparisons.
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(10,5))

ax1.plot(fit_psin, fit_te*lp_te_mult, color=colors[1], lw=lw, zorder=10, label="Input")
ax1.scatter(lp_psin, lp_te, color=colors[2], marker="o", s=s_lp, zorder=11, edgecolor="k", label="LP")
ax1.scatter(dts_psin, dts_te, color=colors[4], marker="^", s=s_lp, zorder=12, edgecolor="k", label="DTS")
ax1.set_xlim(xlims)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.grid()
#ax1.set_xlabel("Psin", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax1.set_ylim([0, None])
ax1.legend(fontsize=fontsize-2)
ax1.set_title("Target Data", fontsize=fontsize)
ax1.text(0.02, 0.9, "a)", fontsize=fontsize, transform=ax1.transAxes)

ax4.plot(fit_psin, fit_jsat*lp_jsat_mult, color=colors[1], lw=lw, label="Input", zorder=10)
ax4.scatter(lp_psin, lp_jsat, color=colors[2], label="LP", marker="o", s=s_lp, zorder=11, edgecolor="k")
ax4.scatter(dts_psin, dts_jsat, color=colors[4], label="DTS", marker="^", s=s_lp, zorder=12, edgecolor="k")
ax4.set_xlim(xlims)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)
ax4.grid()
#ax4.legend(fontsize=fontsize)
ax4.set_xlabel(r"$\mathrm{\psi_N}$", fontsize=fontsize)
ax4.set_ylabel(r"$\mathdefault{j_{sat}\ (A/m^2)}$", fontsize=fontsize)
ax4.ticklabel_format(axis="y", style="scientific", scilimits=(4,4))
ax4.set_ylim([0, None])
ax4.set_yticks(jsat_yticks)
#ax4.text(1.17, 60000, "#167247", rotation=90, bbox=dict(facecolor="white", edgecolor="none"))
ax4.text(0.02, 0.9, "b)", fontsize=fontsize, transform=ax4.transAxes)

ax2.scatter(core_ts["RMRSOMP"] * 100, core_ts["Te (eV)"], s=8, alpha=0.4,
  color=colors[2], marker="^")
ax2.plot(core_ts["RMRSOMP"] * 100, core_ts["OSM Te"], lw=lw, color=colors[1])
ax2.set_xlim([0, 7.5])
ax2.set_ylim([0, 80])
ax2.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax2.grid()
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(axis="x", which="both", labelbottom=False)
ax2.tick_params(axis="y", which="both", labelsize=13)
ax2.set_yticks(np.arange(0, 80, 20))
legend_elements = [Line2D([0], [0], color=colors[1], lw=lw, label="OSM"),
                   Line2D([0], [0], color=colors[2], marker="^", markersize=8, lw=0, label="Thomson")]
ax2.legend(handles=legend_elements, fontsize=fontsize-2)
ax2.set_title("Upstream Comparison", fontsize=fontsize)
ax2.text(0.02, 0.9, "c)", fontsize=fontsize, transform=ax2.transAxes)

ax5.scatter(core_ts["RMRSOMP"] * 100, core_ts["ne (m-3)"], s=8, alpha=0.4,
  color=colors[2], marker="^")
ax5.plot(core_ts["RMRSOMP"] * 100, core_ts["OSM ne"], lw=lw, color=colors[1])
ax5.set_xlim([0, 7.5])
ax5.set_xlabel(r"R-$\mathdefault{R_{sep}}$ OMP (cm)", fontsize=fontsize)
ax5.set_ylabel(r"$\mathdefault{n_e\ (m}^{-\mathdefault{3}})$", fontsize=fontsize)
ax5.grid()
ax5.spines["top"].set_visible(False)
ax5.spines["right"].set_visible(False)
ax5.tick_params(axis="both", which="both", labelsize=fontsize)
ax5.text(0.07, 0.9, "d)", fontsize=fontsize, transform=ax5.transAxes)

ax3.plot(osm_s, osm_te, color=colors[1], lw=lw)
#ax3.scatter(ring_core["S (m)"], ring_core["Te (eV)"], marker=".", s=s,
#  color=colors[2], alpha=0.3, zorder=50)
ax3.scatter(ring_core["S (m)"], ring_core["Te (eV)"], marker="^", s=s_lp,
  color=colors[2], alpha=0.3, zorder=50)
ax3.errorbar(core_avg_s, core_avg_te, yerr=core_std_te, marker="^", ms=8,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[2], zorder=51)
ax3.scatter(ring_div["S (m)"], ring_div["Te (eV)"], marker="^", s=s_lp,
  facecolors=colors[4], edgecolors=colors[4], alpha=0.3, zorder=52)
ax3.errorbar(div_avg_s, div_avg_te, yerr=div_std_te, marker="^", ms=8,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[4], zorder=53)
ax3.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax3.grid()
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.tick_params(axis="x", which="both", labelbottom=False)
ax3.tick_params(axis="y", which="both", labelsize=13)
legend_elements = [Line2D([0], [0], color=colors[1], lw=lw, label="OSM"),
                   Line2D([0], [0], color=colors[2], marker="^", markersize=8, lw=0, label="Core"),
                   #Line2D([0], [0], color=colors[4], marker="^", markersize=8, mfc="none", lw=0, label="Divertor")
                   Line2D([0], [0], color=colors[4], marker="^", markersize=8, lw=0, label="Divertor")]
ax3.legend(handles=legend_elements, fontsize=fontsize-2, loc="upper left")
ax3.set_ylim(ts_te_ylim)
ax3.set_title("Along field line comparison", fontsize=fontsize)
ax3.text(0.5, 0.74, mid_str+"\n"+psin_str, fontsize=fontsize, bbox=dict(color="white"), transform=ax3.transAxes)
ax3.text(0.02, 0.47, "e)", fontsize=fontsize, transform=ax3.transAxes)

ax6.scatter(ring_core["S (m)"], ring_core["ne (m-3)"], marker="^", s=s_lp,
  color=colors[2], alpha=0.3, zorder=50)
ax6.errorbar(core_avg_s, core_avg_ne, yerr=core_std_ne, marker="^", ms=8,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[2], zorder=51)
ax6.scatter(ring_div["S (m)"], ring_div["ne (m-3)"], marker="^", s=s_lp,
  facecolors=colors[4], edgecolors=colors[4], alpha=0.3, zorder=52)
ax6.errorbar(div_avg_s, div_avg_ne, yerr=div_std_ne, marker="^", ms=8,
  color="k", capsize=5, elinewidth=2, mec="k", mfc=colors[4], zorder=53)
ax6.plot(osm_s, osm_ne, color=colors[1], lw=lw)
ax6.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax6.set_ylabel(r"$\mathdefault{n_e\ (m}^{-\mathdefault{3}})$", fontsize=fontsize)
ax6.grid()
ax6.spines["top"].set_visible(False)
ax6.spines["right"].set_visible(False)
ax6.tick_params(axis="both", which="both", labelsize=fontsize)
ax6.text(0.5, 0.74, mid_str+"\n"+psin_str, fontsize=fontsize, bbox=dict(color="white"), transform=ax6.transAxes)
ax6.text(0.02, 0.9, "f)", fontsize=fontsize, transform=ax6.transAxes)

fig.tight_layout()
fig.show()
