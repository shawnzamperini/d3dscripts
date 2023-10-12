# Compare radial profiles of the plasma background to RCP data for 167196.
import pandas as pd
import oedge_plots
import matplotlib.pyplot as plt
import netCDF4
import pickle
import numpy as np


#ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d-167196-modE-shelf-bg.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d-167196-modE-shelf-bg-shifted-smooth-2.nc"
#ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-bg-shifted-ring-entry-06.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-bg-shifted-ring-entry-09.nc"
print("Loading: {}...".format(ncpath))
op = oedge_plots.OedgePlots(ncpath)
print("Calculating fake_probe, could take a while...")
numlocs = 20
rings_only = False  # Turn on for highest resolution, takes a while to run...
print("RCP ne")
op_ne_rcp = op.fake_probe(2.12, 2.37, -0.188, -0.188, data="ne", show_plot=False, plot="psin", num_locs=numlocs, rings_only=rings_only)
print("RCP Te")
op_te_rcp = op.fake_probe(2.20, 2.37, -0.188, -0.188, data="Te", show_plot=False, plot="psin", num_locs=numlocs, rings_only=rings_only)
print("TS ne")
op_ne_ts = op.fake_probe(1.94, 1.94,  0.60,   0.90,  data="ne", show_plot=False, plot="psin", num_locs=numlocs, rings_only=rings_only)
print("TS Te")
op_te_ts = op.fake_probe(1.94, 1.94,  0.60,   0.90,  data="Te", show_plot=False, plot="psin", num_locs=numlocs, rings_only=rings_only)

# Mask out zeros.
#mask = op_ne != 0
#op_psin = op_psin[mask]
#op_ne = op_ne[mask]
#op_te = op_te[mask]

# Load RCP data, bin it.
datapath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/data_for_167196.xlsx"
rcp = pd.read_excel(datapath, sheet_name="rcp_data")
rcp_r = rcp["R (cm)"].values / 100
rcp_psin = rcp["psin"]
rcp_ne = rcp["ne (1e18 m-3)"] * 1e18
rcp_te = rcp["Te (eV)"]

# Use the grid to get an interpolation function of psin(R) at the RCP location.
# This can be commented out, it is only needed once so that we can copy/paste
# the results into the Excel sheet.
#op_psin, op_ne = op.fake_probe(2.15, 2.37, -0.188, -0.188, data="ne", show_plot=False, plot="psin")

# Load TS data from OMFIT.
tspath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/OMFITprofiles_167196_FIT.nc"
ts = netCDF4.Dataset(tspath)
ts_psin = ts["psi_n"][:]
ts_ne = ts["n_e"][:].mean(axis=0)
ts_te = ts["T_e"][:].mean(axis=0)

# I pickled the DTS data with the following in OMFITprofiles:
"""
import pickle

output = {}
dts = OMFIT["OMFITprofiles"]["OUTPUTS"]["RAW"]["TS"]["divertor_r-1"]
for chan in range(1, 8):
    r = dts[chan]["R"]
    z = dts[chan]["Z"]
    psin = dts[chan]["psi_n"]
    te = dts[chan]["T_e"]
    ne = dts[chan]["n_e"]
    output[chan] = {"r":r, "z":z, "psin":psin, "te":te, "ne":ne}
with open("/home/zamperinis/167195_dts.pickle", "wb") as f:
    pickle.dump(output, f)
"""
dtspath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/167195_dts.pickle"
with open(dtspath, "rb") as f:
    udts = pickle.load(f)

# Convert from uncertainty objects to normal values without them.
dts = {}
for chan, d in udts.items():
    te = np.array([v.n for v in udts[chan]["te"][0]])
    ne = np.array([v.n for v in udts[chan]["ne"][0]])
    dts[chan] = {"psin":udts[chan]["psin"], "r":udts[chan]["r"],
        "z":udts[chan]["z"], "te":te, "ne":ne}

# Start of the first window ring.
psin_window = 1.11671590805054

fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, figsize=(9,8), sharex=True)

# RCP comparisons.
ax1.scatter(rcp_psin, rcp_ne, s=15, color="k", alpha=0.75)
ax1.plot(op_ne_rcp["psin"], op_ne_rcp["ne"], color="r")
if rings_only:
    ax1.scatter(op_ne_rcp["psin"], op_ne_rcp["ne"], color="r")
    for i in range(0, len(op_ne_rcp["psin"])):
        ax1.annotate(op_ne_rcp["ring"][i], (op_ne_rcp["psin"][i], op_ne_rcp["ne"][i]))
fig.supxlabel("Psin")
ax1.set_ylabel("ne (m-3)")
ax2.scatter(rcp_psin, rcp_te, s=15, color="k", alpha=0.75, label="RCP")
ax2.plot(op_te_rcp["psin"], op_te_rcp["Te"], color="r", label="OSM")
if rings_only:
    ax2.scatter(op_te_rcp["psin"], op_te_rcp["Te"], color="r")
    for i in range(0, len(op_te_rcp["psin"])):
        ax2.annotate(op_te_rcp["ring"][i], (op_te_rcp["psin"][i], op_te_rcp["Te"][i]))
ax2.set_ylabel("Te (eV)")
ax2.legend()
ax1.set_xlim([0.8, 1.36])
ax1.set_ylim([0, 3e19])
ax2.set_ylim([0, 150])
ax1.axvline(psin_window, color="k", linestyle="--")
ax1.axvline(1.0, color="k", linestyle="-")
ax2.axvline(psin_window, color="k", linestyle="--")
ax2.axvline(1.0, color="k", linestyle="-")
ax1.set_title("RCP")

# TS comparisons.
ax3.plot(ts_psin, ts_ne, color="k", label="TS")
ax3.plot(op_ne_ts["psin"], op_ne_ts["ne"], color="r", label="OSM")
if rings_only:
    ax3.scatter(op_ne_ts["psin"], op_ne_ts["ne"], color="r")
    for i in range(0, len(op_ne_ts["psin"])):
        ax3.annotate(op_ne_ts["ring"][i], (op_ne_ts["psin"][i], op_ne_ts["ne"][i]))
ax4.plot(ts_psin, ts_te, color="k", label="TS")
ax4.plot(op_te_ts["psin"], op_te_ts["Te"], color="r", label="OSM")
if rings_only:
    ax4.scatter(op_te_ts["psin"], op_te_ts["Te"], color="r")
    for i in range(0, len(op_te_ts["psin"])):
        ax4.annotate(op_te_ts["ring"][i], (op_te_ts["psin"][i], op_te_ts["Te"][i]))
ax3.set_ylabel("ne (m-3)")
ax4.set_ylabel("Te (eV)")
ax4.legend()
ax3.set_ylim([0, 3e19])
ax4.set_ylim([0, 150])
ax3.axvline(psin_window, color="k", linestyle="--")
ax3.axvline(1.0, color="k", linestyle="-")
ax4.axvline(psin_window, color="k", linestyle="--")
ax4.axvline(1.0, color="k", linestyle="-")
ax3.set_title("Thomson Scattering")

fig.tight_layout()
fig.show()
