# Compare a SOL29 run to experimentally measured data.
import get_lp
from ThomsonClass import ThomsonClass
import matplotlib.pyplot as plt
import numpy as np
import oedge_plots
import pandas as pd

import sys
sys.path.append("/usr/local/mdsplus/python")


# The shot we are looking at.
shot = 167195
time = 2500
tmin = 3800
tmax = 5000
tmany = 15
tstart = tmin
tend = tmax

# Load LP data.
lpdict = get_lp.plot_lps(shot, tstart, tend, bins=100, tunnel=False, showplot=False)
labels = ["S-1", "S-2", "S-3", "S-4", "S-5", "S-6", "S-7", "S-8", "S-9",
    "S10", "S11"]

colors = {}
for i in range(0, len(labels)):
    colors[labels[i]] = "C{}".format(i)

data = {"rmrs":[], "ne":[], "te":[], "color":[]}
for i in range(0, len(lpdict["labels"])):

    label = lpdict["labels"][i].strip()
    if label in labels:
        data["rmrs"].append(lpdict["rminrsep"][i])
        data["ne"].append(lpdict["ne (cm-3)"][i] * 1e6)
        data["te"].append(lpdict["Te (eV)"][i])
        data["color"].append(colors[label])

sort_idx = np.argsort(data["rmrs"])
data["rmrs"] = np.array(data["rmrs"])[sort_idx]
data["ne"] = np.array(data["ne"])[sort_idx]
data["te"] = np.array(data["te"])[sort_idx]
data["color"] = np.array(data["color"])[sort_idx]

# Load the TS data.
ts = ThomsonClass(shot, 'core')
ts.load_ts(tunnel=False)
ts.map_to_efit(np.linspace(tmin, tmax, tmany), tree="EFIT01")

# Pull out the arrays.
r  = ts.avg_omp['RminRsep_omp'] * 100
te = ts.avg_omp['Te_omp']
ne = ts.avg_omp['ne_omp'] * 10**(-18)
r_err  = ts.avg_omp['RminRsep_omp_err'] * 100  # m to cm
te_err = ts.avg_omp['Te_omp_err']
ne_err = ts.avg_omp['ne_omp_err'] * 10**(-18) # m-3 to 10^-18 m-3

# The Z values and average value between time range of interest.
ts_z = ts.ts_dict["z"]["Y"]
ts_time = ts.ts_dict["temp"]["X"]
mask = np.logical_and(ts_time>=tmin, ts_time<=tmax)
ts_te = ts.ts_dict["temp"]["Y"][:,mask]
ts_ne = ts.ts_dict["density"]["Y"][:,mask]

# Load the RCP data.
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/MP167195_2.tab"
rcp = pd.read_csv(rcp_path, delimiter="\t")
rcp_r = rcp["R(cm)"].values
rcp_te = rcp["Te(eV)"].values
rcp_ne = rcp["Ne(E18 m-3)"].values * 1e18

# Load the SOL29 results at each location we care about.
ncpath  = "/Users/zamperini/Documents/d3d_work/divimp_files/sol29_testing/d3d-sol29-test-001.nc"
op = oedge_plots.OedgePlots(ncpath)

# Load target data, restrict to just the second target (the outer here).
outer_end = int(op.nc["NDSIN"][:])
op_targ_x = op.nc["SEPDIS"][:][:outer_end]
op_targ_te = op.nc["KTEDS"][:][:outer_end]
op_targ_ne = op.nc["KNDS"][:][:outer_end]
op_targ_counts = op.nc["blob_counts_targ"][:][:outer_end]

nks = np.array(op.nc["NKS"[:]]) - 1
op_targ_x = []
op_targ_te = []
op_targ_ne = []
op_targ_counts = []
for ring in range(18, 70):
    op_targ_x.append(float(op.nc["RS"][ring, nks[ring]]))
    op_targ_te.append(float(op.nc["KTEBS"][ring, nks[ring]]))
    op_targ_ne.append(float(op.nc["KNBS"][ring, nks[ring]]))
    op_targ_counts.append(float(op.nc["blob_counts"][ring, nks[ring]]))

# Convert to R-Rsep.
rsep = 1.426
op_targ_x = np.array(op_targ_x) - rsep

# Pull out a radial profile where the RCP plunges.
op_tsc_te = op.fake_probe(1.94, 1.94, 0.70, 0.85, data="Te", plot="psin", show_plot=False)
op_tsc_ne = op.fake_probe(1.94, 1.94, 0.70, 0.85, data="ne", plot="psin", show_plot=False)
op_tsd_te = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="Te", plot="psin", show_plot=False)
op_tsd_ne = op.fake_probe(1.484, 1.484, -0.82, -1.17, data="ne", plot="psin", show_plot=False)
op_rcp_te = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Te", plot="R", show_plot=False)
op_rcp_ne = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="ne", plot="R", show_plot=False)
op_rcp_m  = op.fake_probe(2.18, 2.30, -0.188, -0.188, data="Mach", plot="R", show_plot=False)


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(10,6))

ax1.axvline(0, color="k", linestyle="--")
ax1.scatter(data["rmrs"], data["te"], c=data["color"], zorder=5)
ax1.plot(op_targ_x, op_targ_te, zorder=10)
ax1.set_xlabel("Distance from separatrix")
ax1.set_ylabel("Te (eV)")
ax1.set_xlim([-0.05, 0.25])

ax44 = ax4.twinx()
ax4.axvline(0, color="k", linestyle="--")
ax4.scatter(data["rmrs"], data["ne"], c=data["color"], zorder=5)
ax4.plot(op_targ_x, op_targ_ne, zorder=10)
ax44.plot(op_targ_x, op_targ_counts, zorder=15, color="r")
ax4.set_xlabel("Distance from separatrix")
ax4.set_ylabel("LP ne (m-3)")
ax44.set_ylabel("Blob Counts", color="r")
ax4.set_xlim([-0.05, 0.25])

ax2.axvline(2.2169, color="k", linestyle="--")
ax2.scatter(rcp_r/100, rcp_te)
ax2.plot(op_rcp_te["r"], op_rcp_te["Te"])

ax5.axvline(2.2169, color="k", linestyle="--")
ax5.scatter(rcp_r/100, rcp_ne)
ax5.plot(op_rcp_ne["r"], op_rcp_ne["ne"])
ax5.set_ylim([0, 3e19])

ax3.axvline(0.733, color="k", linestyle="--")
ax6.axvline(0.733, color="k", linestyle="--")
ts_z_all = np.repeat(ts_z, ts_te.shape[1])
ts_te_all = ts_te.flatten()
ts_ne_all = ts_ne.flatten()
#for i in range(0, ts_te.shape[1]):
ax3.scatter(ts_z_all, ts_te_all, color="k", alpha=0.2, marker=".")
ax6.scatter(ts_z_all, ts_ne_all, color="k", alpha=0.2, marker=".")

ax3.plot(op_tsc_te["z"], op_tsc_te["Te"])
ax3.set_xlim([0.7, 0.85])
ax3.set_ylim([0, 150])

ax6.plot(op_tsc_ne["z"], op_tsc_ne["ne"])
ax6.set_xlim([0.7, 0.85])
ax6.set_ylim([0, 2.5e19])

fig.tight_layout()
fig.show()
