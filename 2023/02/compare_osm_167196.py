# Script to compare the background of 167196 to some of the available diagnostics
# from that day.
import oedge_plots
import matplotlib.pyplot as plt
import netCDF4
import pickle
import numpy as np
import pandas as pd
from scipy.signal import medfilt
from scipy.interpolate import griddata
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

sys.path.append("../../2022/12/")
import BlobbyFarSOL

# Load background.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-bg-shifted-ring-entry-13.nc"
ncpath_mach = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-017.nc"
op = oedge_plots.OedgePlots(ncpath)
op_mach = oedge_plots.OedgePlots(ncpath_mach)

# To copy/paste into OMFIT for getting the below data with the fastTS module.
"""
import pickle
import numpy as np

root = OMFIT['fastTS']['OUTPUTS']['current']['filtered_TS']
shot = int(OMFIT['fastTS']['OUTPUTS']['current']['filtered_TS']['shot'])

output = {}
for sysname in ["core", "divertor", "tangential"]:
    sys = root[sysname]
    tmp = {}
    tmp["time"] = np.array(sys["time"])
    tmp["r"] = np.array(sys["r"])
    tmp["z"] = np.array(sys["z"])
    tmp["te"] = np.array(sys["temp"])
    tmp["ne"] = np.array(sys["density"])
    tmp["te_err"] = np.array(sys["temp_e"])
    tmp["ne_err"] = np.array(sys["density_e"])
    tmp["psin"] = np.array(sys["psin_TS"])
    tmp["chord"] = np.array(sys["chord_index"])

    output[sysname] = tmp

with open("/home/zamperinis/ts_{}.pickle".format(shot), "wb") as f:
    pickle.dump(output, f)

"""

# Load Thomson scattering data from OMFIT. Put into ts_plot just the data during the flattop.
corets_shift = 0.03
print("Core TS shifted by psin of {}".format(corets_shift))
ts_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/ts_167196.pickle"
with open(ts_path, "rb") as f:
    ts = pickle.load(f)
ts_plot = {"core":{}, "divertor":{}, "tangential":{}}
for sys in ts.keys():
    tmp = ts[sys]
    mask = np.logical_and(tmp["time"]>=2500, tmp["time"]<=5000)
    ts_plot[sys]["time"] = tmp["time"][mask]
    for key in ["te", "te_err", "ne", "ne_err", "psin"]:
        if sys == "core" and key == "psin":
            ts_plot[sys][key] = tmp[key][:, mask] + corets_shift
        else:
            ts_plot[sys][key] = tmp[key][:,mask]
    ts_plot[sys]["chord"] = tmp["chord"]

# Load LLAMA data.
# llama = np.load("/Users/zamperini/Documents/d3d_work/files/LLAMA_190484_.npz")
# lpsin = llama["psi_n"]
# lneut = medfilt(llama["nDens_LFS"], 25)
# lneut_err = llama["nDens_LFS_err"]
# lion = llama["ion_LFS"]
# lion_err = llama["ion_LFS_err"]

# Get some gfile stuff so we can go from R, Z to psin.
gfile_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/167195_4460.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
psin_2D = gfile["PSIRZ_NORM"]

# Load RCP data. Map to psin.
rcp_shift = -0.015
print("RCP shifted by {:.3} m".format(rcp_shift))
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP167195_2.tab"
rcp = pd.read_csv(rcp_path, delimiter="\t")
# rcp_coord = zip((rcp["R(cm)"].values + rcp_shift) / 100, np.full(len(rcp["R(cm)"].values), -0.185))
rcp_coord = (rcp["R(cm)"].to_numpy() / 100 + rcp_shift, np.full(len(rcp["R(cm)"].values), -0.185))
rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), rcp_coord)

# Do fake plunges at each diagnostic location for OSM comparisons.
# TS (core): R = 1.94, Z = 0.50 - 0.72
# TS (div):  R = 1.484, Z = -0.82 - -1.17
# RCP: R = 2.18 - 2.30, Z = -0.188
op_tsc_te = op.fake_probe(1.94, 1.94, 0.67, 0.85, data="Te", plot="psin", show_plot=False)
op_tsc_ne = op.fake_probe(1.94, 1.94, 0.67, 0.85, data="ne", plot="psin", show_plot=False)
op_tsd_te = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="Te", plot="psin", show_plot=False)
op_tsd_ne = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="ne", plot="psin", show_plot=False, rings_only=False)
op_rcp_te = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Te", plot="psin", show_plot=False, rings_only=False)
op_rcp_ne = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="ne", plot="psin", show_plot=False, rings_only=False)
op_rcp_m  = op_mach.fake_probe(2.18, 2.30, -0.188, -0.188, data="Mach", plot="psin", show_plot=False, rings_only=True)
div_llama = op.fake_probe(1.93, 2.05, -0.77, -0.77, data="neut_dens", plot="psin", show_plot=False)
try:
    div_rcp   = op.fake_probe(2.23, 2.33, -0.185, -0.185, data="neut_dens", plot="psin", show_plot=False)
    neut = True
except:
    neut = False

# Can include neutral estimates from the blob analysis at least.
if neut:
    bfs = BlobbyFarSOL.main(167195, 2, showplot=False, temod=0.50)
    blob_coord = (bfs.blob_r / 100 + rcp_shift, np.full(len(bfs.blob_r), -0.185))
    blob_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), blob_coord)

# Pull out ring, psin values so we can plot vertical lines at them.
ring_linex = []; ring_liney = []
# for i in range(0, len(op_tsc_te["ring"])):
#     if op_tsc_te["ring"][i] in [10, 20, 30, 40]:
#         ring_linex.append(op_tsc_te["psin"][i])
#         ring_liney.append(op_tsc_te["ring"][i])
for i in range(0, len(op.nc["PSIFL"][:, 0])):
    ring_linex.append(op.nc["PSIFL"][i, 0])
    ring_liney.append(i + 1)

# Now we do our comparison plots.
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(10, 5))

# Vertical lines to orient wrt ring numbers.
for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]:
    for psin, ring in zip(ring_linex, ring_liney):
        if ring in [20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:  # Only the SOL rings.
            ax.axvline(psin, color="k", linestyle="--")
    ax.axvline(1.0, color="k")

# Core TS Te.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["te"].flatten()
yerr = ts_plot["core"]["te_err"].flatten()
ax1.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax1.plot(op_tsc_te["psin"], op_tsc_te["Te"], color="tab:red")
ax1.set_xlabel("Psin")
ax1.set_title("Core TS Te")
ax1.set_xlim([0.99, 1.15])
ax1.set_ylim([0, 100])

# Core TS ne.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["ne"].flatten()
yerr = ts_plot["core"]["ne_err"].flatten()
ax2.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax2.plot(op_tsc_ne["psin"], op_tsc_ne["ne"], color="tab:red")
ax2.set_xlabel("Psin")
ax2.set_title("Core TS ne")
ax2.set_xlim([0.99, 1.15])
ax2.set_ylim([0, 2.0e19])

# Divertor TS Te at the crown (chords 9-13).
# mask = np.logical_and(ts_plot["divertor"]["chord"]>=9, ts_plot["divertor"]["chord"]<=13)
x = ts_plot["divertor"]["psin"].flatten()
y = ts_plot["divertor"]["te"].flatten()
yerr = ts_plot["divertor"]["te_err"].flatten()
ax3.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax3.plot(op_tsd_te["psin"], op_tsd_te["Te"], color="tab:red")
ax3.set_xlabel("Psin")
ax3.set_title("Divertor TS Te")
ax3.set_xlim([0.99, 1.03])
ax3.set_ylim([0, 100])

# Divertor TS ne at the crown.
x = ts_plot["divertor"]["psin"].flatten()
y = ts_plot["divertor"]["ne"].flatten()
yerr = ts_plot["divertor"]["ne_err"].flatten()
ax4.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax4.plot(op_tsd_ne["psin"], op_tsd_ne["ne"], color="tab:red")
ax4.set_xlabel("Psin")
ax4.set_title("Divertor TS ne")
ax4.set_xlim([0.99, 1.03])
ax4.set_ylim([0, 1e20])

# RCP Te.
x = rcp_psin
y = rcp["Te(eV)"].values
ax5.scatter(x, y, s=15, color="k")
ax5.plot(op_rcp_te["psin"], op_rcp_te["Te"], color="tab:red", marker=".")
ax5.set_xlabel("Psin")
ax5.set_title("RCP Te")
# ax5.axvline(2.2367, color="k", linestyle="--")
ax5.set_xlim([0.99, 1.3])
ax5.set_ylim([0, 50])

# RCP ne.
x = rcp_psin
y = rcp["Ne(E18 m-3)"].values * 1e18
ax6.scatter(x, y, s=15, color="k")
ax6.plot(op_rcp_ne["psin"], op_rcp_ne["ne"], color="tab:red", marker=".")
ax6.set_xlabel("Psin")
ax6.set_title("RCP ne")
# ax6.axvline(2.2367, color="k", linestyle="--")
ax6.set_xlim([0.99, 1.3])
ax6.set_ylim([0, 2e19])

# RCP Mach.
x = rcp_psin
y = rcp["Machn"].values
ax7.axhline(0, color="k")
ax7.scatter(x, y, s=15, color="k")
ax7.plot(op_rcp_m["psin"], op_rcp_m["Mach"], color="tab:red")
ax7.set_xlabel("Psin")
ax7.set_title("RCP Mach")
# ax7.axvline(2.2367, color="k", linestyle="--")
ax7.set_xlim([0.99, 1.25])

if neut:
    # Neutral density.
    # ax8.plot(lpsin, lneut, color="k")
    # ax8.plot(div_llama["psin"], div_llama["neut_dens"], color="tab:red")
    ax8.scatter(blob_psin, bfs.neut_dens2, marker="*", s=75, color="tab:cyan", edgecolors="k")
    ax8.plot(div_rcp["psin"], div_rcp["neut_dens"], color="tab:cyan")
    ax8.set_xlabel("Psin")
    ax8.set_title("Neutral Density")
    ax8.set_xlim([0.99, 1.25])
    ax8.set_ylim([1e14, 1e17])
    ax8.set_yscale("log")
    ax8.grid(alpha=0.3)

fig.tight_layout()
fig.show()

# Additional plot comparing the TS to the RCP.
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
#
# ax1.axvline(1.0, color="r")
# ax2.axvline(1.0, color="r")
# x = ts_plot["core"]["psin"].flatten()
# y = ts_plot["core"]["te"].flatten()
# yerr = ts_plot["core"]["te_err"].flatten()
# ax1.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
# x = rcp_psin
# y = rcp["Te(eV)"].values
# ax1.scatter(x, y, s=15, color="k")
# ax1.set_ylim([0, 200])
# ax1.set_ylabel("Te (eV)")
#
# x = ts_plot["core"]["psin"].flatten()
# y = ts_plot["core"]["ne"].flatten()
# yerr = ts_plot["core"]["ne_err"].flatten()
# ax2.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
# x = rcp_psin
# y = rcp["Ne(E18 m-3)"].values * 1e18
# ax2.scatter(x, y, s=15, color="k")
# ax2.set_xlim([0.95, 1.3])
# ax2.set_ylim([0, 2.5e19])
# ax2.set_ylabel("ne (1e18 m-3)")
#
# fig.tight_layout()
# fig.show()

# An additional pared down plot for the paper.
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(7, 5))

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.axvline(1.0, color="k", linestyle="--")

# Core TS Te.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["te"].flatten()
yerr = ts_plot["core"]["te_err"].flatten()
ax1.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0, alpha=0.7)
ax1.plot(op_tsc_te["psin"], op_tsc_te["Te"], color="tab:pink", label=r"TS $\mathdefault{T_e}$ (eV)")
ax1.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax1.set_title("Core TS Te")
ax1.set_xlim([0.97, 1.15])
ax1.set_ylim([0, 120])

# Core TS ne.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["ne"].flatten()
yerr = ts_plot["core"]["ne_err"].flatten()
ax2.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0, alpha=0.7)
ax2.plot(op_tsc_ne["psin"], op_tsc_ne["ne"], color="tab:pink", label=r"TS $\mathdefault{n_e\ (m^{-3})}$")
ax2.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax2.set_title("Core TS ne")
ax2.set_xlim([0.97, 1.15])
ax2.set_ylim([0, 2.5e19])
ax2.set_yticks([0, 1e19, 2e19])

# RCP Te.
x = rcp_psin
y = rcp["Te(eV)"].values
ax4.scatter(x, y, s=15, color="k")
ax4.plot(op_rcp_te["psin"], op_rcp_te["Te"], color="tab:pink", label=r"RCP $\mathdefault{T_e}$ (eV)")
ax4.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax4.set_title("RCP Te")
# ax4.axvline(2.2367, color="k", linestyle="--")
ax4.set_xlim([0.99, 1.3])
ax4.set_ylim([0, 70])

# RCP ne.
x = rcp_psin
y = rcp["Ne(E18 m-3)"].values * 1e18
ax5.scatter(x, y, s=15, color="k")
ax5.plot(op_rcp_ne["psin"], op_rcp_ne["ne"], color="tab:pink", label=r"RCP $\mathdefault{n_e\ (m^{-3})}$")
ax5.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax5.set_title("RCP ne")
# ax5.axvline(2.2367, color="k", linestyle="--")
ax5.set_xlim([0.99, 1.3])
ax5.set_ylim([0, 2.5e19])
ax5.set_yticks([0, 1e19, 2e19])

# RCP Mach.
x = rcp_psin
y = rcp["Machn"].values * -1
ax6.axhline(0, color="k")
ax6.scatter(x, y, s=15, color="k")
ax6.plot(op_rcp_m["psin"], op_rcp_m["Mach"], color="tab:pink", label="RCP Mach")
ax6.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax6.set_title("RCP Mach")
# ax6.axvline(2.2367, color="k", linestyle="--")
ax6.set_xlim([0.99, 1.3])
ax6.set_ylim([-1, 0.5])
ax6.set_yticks([-1, -0.5, 0, 0.5])

# Neutral density.
if neut:
    # ax3.scatter(blob_psin, bfs.neut_dens2, marker="*", s=75, color="tab:pink", edgecolors="k")
    ax3.plot(div_rcp["psin"], div_rcp["neut_dens"], color="tab:pink", label=r"OMP $\mathdefault{n_n\ (m^{-3})}$")
    ax3.set_xlabel(r"$\mathdefault{\psi_n}$")
    # ax3.set_title("Neutral Density")
    ax3.set_xlim([0.99, 1.3])
    ax3.set_ylim([1e15, 1e18])
    ax3.set_yscale("log")
    # ax3.grid(alpha=0.3)

axs = [ax1, ax2, ax3, ax4, ax5, ax6]
letters = ["a", "b", "c", "d", "e", "f"]
for i in range(0, len(axs)):
    axs[i].legend(fontsize=10, loc="upper right")
    axs[i].grid(alpha=0.3)
    axs[i].text(0.04, 0.04, "{})".format(letters[i]), transform=axs[i].transAxes)

fig.tight_layout()
fig.show()

# A smaller plot that I can use in my 4-pager SiC paper.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3, 4.5))
ax1.axvline(1.0, color="k", linestyle="--")
ax2.axvline(1.0, color="k", linestyle="--")

rcp_mask = rcp_psin < 1.24
x = rcp_psin[rcp_mask]
y = rcp["Ne(E18 m-3)"].values[rcp_mask] * 1e18
ax1.scatter(x, y, s=12, color="k")
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["ne"].flatten()
yerr = ts_plot["core"]["ne_err"].flatten()
ax1.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0, alpha=0.7)
ax1.plot(op_rcp_ne["psin"], op_rcp_ne["ne"], color="k", lw=3)
ax1.plot(op_rcp_ne["psin"], op_rcp_ne["ne"], color="tab:cyan", label=r"RCP $\mathdefault{n_e\ (m^{-3})}$", lw=2)
# ax1.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax5.set_title("RCP ne")
# ax5.axvline(2.2367, color="k", linestyle="--")
ax1.set_xlim([0.9, 1.3])
ax1.set_ylim([1.0e18, 4.0e19])
#ax1.set_yticks([0, 1e19, 2e19])
#ax1.legend()
ax1.set_yscale("log")
ax1.set_ylabel(r"$\mathdefault{n_e\ (m^{-3})}$")
ax1.text(1.01, 0.7, "#167196", transform=ax1.transAxes, rotation=270, fontsize=8)

# Inset of the plasma shape.
axins = inset_axes(ax1, width=0.6, height=0.8)
axins.contour(Rs, Zs, psin_2D, levels=[1.0], colors="k", linewidths=1)
axins.set_aspect("equal")
axins.plot(gfile["RLIM"], gfile["ZLIM"], color="k", linewidth=1)
# axins.spines["top"].set_visible(False)
# axins.spines["right"].set_visible(False)
# axins.spines["left"].set_visible(False)
# axins.spines["bottom"].set_visible(False)
axins.tick_params(axis="both", which="both", labelleft=False, labelbottom=False, left=False, bottom=False)

x = rcp_psin[rcp_mask]
y = rcp["Te(eV)"].values[rcp_mask]
ax2.scatter(x, y, s=12, color="k")
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["te"].flatten()
yerr = ts_plot["core"]["te_err"].flatten()
ax2.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0, alpha=0.7)
ax2.plot(op_rcp_te["psin"], op_rcp_te["Te"], color="k", lw=3)
ax2.plot(op_rcp_te["psin"], op_rcp_te["Te"], color="tab:red", label=r"RCP $\mathdefault{T_e}$ (eV)", lw=2)
ax2.set_xlabel(r"$\mathdefault{\psi_n}$")
# ax4.set_title("RCP Te")
# ax4.axvline(2.2367, color="k", linestyle="--")
ax2.set_xlim([0.9, 1.3])
ax2.set_ylim([1, 500])
#ax2.legend()
ax2.set_yscale("log")
ax2.set_ylabel(r"$\mathdefault{T_e\ (eV)}$")

fig.tight_layout()
fig.show()