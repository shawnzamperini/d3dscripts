# Script to compare the background of 190423 to some of the available diagnostics
# from that day.
import oedge_plots
import matplotlib.pyplot as plt
import netCDF4
import pickle
import numpy as np
import pandas as pd


# Load background.
#ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-bkg-002.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-sput-004.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-tungsten-001.nc"
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
    tmp["density_e"] = np.array(sys["density_e"])
    tmp["psin"] = np.array(sys["psin_TS"])
    tmp["chord"] = np.array(sys["chord_index"])

    output[sysname] = tmp

with open("/home/zamperinis/ts_{}.pickle".format(shot), "wb") as f:
    pickle.dump(output, f)

"""

# Load Thomson scattering data from OMFIT. 190422-423 are repeats. Put into
# ts_plot just the data during the flattop.
ts_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/ts_190423.pickle"
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
# To-do.

# Load RCP data. Diagnostic shot technically would be 190419, but forgot to plunge.
# We don't technically have matching RCP data since it was all plunged during
# density ramps. The closest thing we have to a diagnostic plunge would be
# MP190411_2. During the plunge it was ~similar to 190423. Crown was a smidgen
# different but strike point, power are the same, and density is pretty close.
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190411_2.tab"
rcp = pd.read_csv(rcp_path, delimiter="\t")

# Do fake plunges at each diagnostic location for OSM comparisons.
# TS (core): R = 1.94, Z = 0.50 - 0.72
# TS (div):  R = 1.484, Z = -0.82 - -1.17
# RCP: R = 2.18 - 2.30, Z = -0.188
op_tsc_te = op.fake_probe(1.94, 1.94, 0.50, 0.72, data="Te", plot="psin", show_plot=False)
op_tsc_ne = op.fake_probe(1.94, 1.94, 0.50, 0.72, data="ne", plot="psin", show_plot=False)
op_tsd_te = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="Te", plot="psin", show_plot=False)
op_tsd_ne = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="ne", plot="psin", show_plot=False)
op_rcp_te = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Te", plot="R", show_plot=False)
op_rcp_ne = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="ne", plot="R", show_plot=False)
op_rcp_m  = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Mach", plot="R", show_plot=False)


# Pull out ring, psin values so we can plot vertical lines at them.
ring_linex = []; ring_liney = []
for i in range(0, len(op_tsc_te["ring"])):
    if op_tsc_te["ring"][i] in [10, 20, 30, 40]:
        ring_linex.append(op_tsc_te["psin"][i])
        ring_liney.append(op_tsc_te["ring"][i])


# Now we do our comparison plots.
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(10, 5))

# Vertical lines to oreint wrt ring numbers.
for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    for psin, ring in zip(ring_linex, ring_liney):
        ax.axvline(psin, color="k", linestyle="--")

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
ax2.set_ylim([0, 2.5e19])

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
x = rcp["R(cm)"].values / 100
y = rcp["Te(eV)"].values
ax5.scatter(x, y, s=15, color="k")
ax5.plot(op_rcp_te["r"], op_rcp_te["Te"], color="tab:red")
ax5.set_xlabel("R (m)")
ax5.set_title("RCP Te")
ax5.axvline(2.2367, color="k", linestyle="--")
ax5.set_xlim([2.23, 2.36])
ax5.set_ylim([0, 50])

# RCP ne.
x = rcp["R(cm)"].values / 100
y = rcp["Ne(E18 m-3)"].values * 1e18
ax6.scatter(x, y, s=15, color="k")
ax6.plot(op_rcp_ne["r"], op_rcp_ne["ne"], color="tab:red")
ax6.set_xlabel("R (m)")
ax6.set_title("RCP ne")
ax6.axvline(2.2367, color="k", linestyle="--")
ax6.set_xlim([2.23, 2.36])
#ax5.set_ylim([0, 50])

# RCP Mach.
x = rcp["R(cm)"].values / 100
y = rcp["Machn"].values
ax7.axhline(0, color="k")
ax7.scatter(x, y, s=15, color="k")
ax7.plot(op_rcp_m["r"], op_rcp_m["Mach"], color="tab:red")
ax7.set_xlabel("R (m)")
ax7.set_title("RCP Mach")
ax7.axvline(2.2367, color="k", linestyle="--")
ax7.set_xlim([2.23, 2.36])

fig.tight_layout()
fig.show()
