# Script to compare the background of 190423 to some of the available diagnostics
# from that day.
import oedge_plots
import matplotlib.pyplot as plt
import netCDF4
import pickle
import numpy as np
import pandas as pd
from scipy.signal import medfilt
from scipy.interpolate import griddata
import h5py


# Load background.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190484/d3d-190484-bkg-005-outgas-smooth.nc"
op = oedge_plots.OedgePlots(ncpath)

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
ts_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190484/ts_190484.pickle"
with open(ts_path, "rb") as f:
    ts = pickle.load(f)
ts_plot = {"core":{}, "divertor":{}, "tangential":{}}
for sys in ts.keys():
    tmp = ts[sys]
    mask = np.logical_and(tmp["time"]>=2500, tmp["time"]<=5000)
    ts_plot[sys]["time"] = tmp["time"][mask]
    for key in ["te", "te_err", "ne", "ne_err", "psin"]:
        ts_plot[sys][key] = tmp[key][:,mask]
    ts_plot[sys]["chord"] = tmp["chord"]

# Load LLAMA data.
llama = np.load("/Users/zamperini/Documents/d3d_work/files/LLAMA_190484_.npz")
lpsin = llama["psi_n"]
lneut = medfilt(llama["nDens_LFS"], 25)
lneut_err = llama["nDens_LFS_err"]
lion = llama["ion_LFS"]
lion_err = llama["ion_LFS_err"]

# Get some gfile stuff so we can go from R, Z to psin.
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]

# Load RCP data. Map to psin.
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190484_1.tab"
rcp = pd.read_csv(rcp_path, delimiter="\t")
rcp_coord = zip(rcp["R(cm)"].values / 100, np.full(len(rcp["R(cm)"].values), -0.185))
rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))

# Load Kirtan's estimate of the neutral density from the filterscopes.
filt_path = "/Users/zamperini/My Drive/Research/Data/Neutral_Density_and_Ionization_Rate_190484_Hmode.h5"
filt = h5py.File(filt_path)
filt_rho = filt["Rho_spline"][:]
filt_neut = filt["Neutral_Density_Da"][:] * 1e6 # cm-3 to m-3
filt_neut_err = filt["Neutral_Density_Da_Err"][:]
filt_psin = np.sqrt(filt_rho)

# Do fake plunges at each diagnostic location for OSM comparisons.
# TS (core): R = 1.94, Z = 0.50 - 0.72
# TS (div):  R = 1.484, Z = -0.82 - -1.17
# RCP: R = 2.18 - 2.30, Z = -0.188
op_tsc_te = op.fake_probe(1.94, 1.94, 0.50, 0.72, data="Te", plot="psin", show_plot=False)
op_tsc_ne = op.fake_probe(1.94, 1.94, 0.50, 0.72, data="ne", plot="psin", show_plot=False)
op_tsd_te = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="Te", plot="psin", show_plot=False)
op_tsd_ne = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="ne", plot="psin", show_plot=False)
op_rcp_te = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Te", plot="psin", show_plot=False, rings_only=True)
op_rcp_ne = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="ne", plot="psin", show_plot=False, rings_only=True)
op_rcp_m  = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Mach", plot="psin", show_plot=False)
div_llama = op.fake_probe(1.93, 2.05, -0.77, -0.77, data="neut_dens", plot="psin", show_plot=False)
div_rcp   = op.fake_probe(2.23, 2.33, -0.185, -0.185, data="neut_dens", plot="psin", show_plot=False)


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
        if ring in [20, 30, 40]:  # Only the SOL rings.
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
ax1.set_xlim([0.99, 1.11])
ax1.set_ylim([0, 100])

# Core TS ne.
x = ts_plot["core"]["psin"].flatten()
y = ts_plot["core"]["ne"].flatten()
yerr = ts_plot["core"]["ne_err"].flatten()
ax2.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax2.plot(op_tsc_ne["psin"], op_tsc_ne["ne"], color="tab:red")
ax2.set_xlabel("Psin")
ax2.set_title("Core TS ne")
ax2.set_xlim([0.99, 1.11])
ax2.set_ylim([0, 2.0e19])

# Divertor TS Te at the crown (chords 9-13).
mask = np.logical_and(ts_plot["divertor"]["chord"]>=9, ts_plot["divertor"]["chord"]<=13)
x = ts_plot["divertor"]["psin"][mask].flatten()
y = ts_plot["divertor"]["te"][mask].flatten()
yerr = ts_plot["divertor"]["te_err"][mask].flatten()
ax3.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax3.plot(op_tsd_te["psin"], op_tsd_te["Te"], color="tab:red")
ax3.set_xlabel("Psin")
ax3.set_title("Divertor TS Te")
ax3.set_xlim([0.99, 1.11])
ax3.set_ylim([0, 100])

# Divertor TS ne at the crown.
x = ts_plot["divertor"]["psin"][mask].flatten()
y = ts_plot["divertor"]["ne"][mask].flatten()
yerr = ts_plot["divertor"]["ne_err"][mask].flatten()
ax4.errorbar(x, y, yerr, elinewidth=1, ecolor="k", color="k", markersize=15, lw=0)
ax4.plot(op_tsd_ne["psin"], op_tsd_ne["ne"], color="tab:red")
ax4.set_xlabel("Psin")
ax4.set_title("Divertor TS ne")
ax4.set_xlim([0.99, 1.11])
ax4.set_ylim([0, 2.5e19])

# RCP Te.
x = rcp_psin
y = rcp["Te(eV)"].values
ax5.scatter(x, y, s=15, color="k")
ax5.plot(op_rcp_te["psin"], op_rcp_te["Te"], color="tab:red", marker=".")
ax5.set_xlabel("Psin")
ax5.set_title("RCP Te")
# ax5.axvline(2.2367, color="k", linestyle="--")
ax5.set_xlim([0.99, 1.25])
ax5.set_ylim([0, 50])

# RCP ne.
x = rcp_psin
y = rcp["Ne(E18 m-3)"].values * 1e18
ax6.scatter(x, y, s=15, color="k")
ax6.plot(op_rcp_ne["psin"], op_rcp_ne["ne"], color="tab:red", marker=".")
ax6.set_xlabel("Psin")
ax6.set_title("RCP ne")
# ax6.axvline(2.2367, color="k", linestyle="--")
ax6.set_xlim([0.99, 1.25])
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

# Neutral density.
ax8.plot(lpsin, lneut, color="k")
ax8.plot(filt_psin, filt_neut, color="r")
ax8.plot(div_llama["psin"], div_llama["neut_dens"], color="tab:red")
ax8.plot(div_rcp["psin"], div_rcp["neut_dens"], color="tab:cyan")
ax8.set_xlabel("Psin")
ax8.set_title("Neutral Density")
ax8.set_xlim([0.99, 1.11])
ax8.set_ylim([1e14, 1e17])
ax8.set_yscale("log")
ax8.grid(alpha=0.3)

fig.tight_layout()
fig.show()
