# Plot results from a scan in correlation time.
import oedge_plots
import LimPlots
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd


root = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div20"
cases = ["a", "b", "c", "d", "e"]
corrs = [5, 10, 15, 20, 25]
ops = []
for case in tqdm(cases):
    path = root + case + ".nc"
    ops.append(oedge_plots.OedgePlots(path))

ss_near = []; nzs_near = []
ss_far = []; nzs_far = []
for op in tqdm(ops):
    s, nz = op.along_ring(91, "DDLIMS", charge="all", plot_it=False)
    ss_far.append(s)
    nzs_far.append(nz)
    s, nz = op.along_ring(22, "DDLIMS", charge="all", plot_it=False)
    ss_near.append(s)
    nzs_near.append(nz)

near_rmrsomp = float(op.nc["MIDIST"][1][22+1])
far_rmrsomp = float(op.nc["MIDIST"][1][91+1])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
for i in range(0, len(ops)):
    ax1.plot(ss_near[i], nzs_near[i], label=corrs[i])
    ax2.plot(ss_far[i], nzs_far[i], label=corrs[i])
ax1.grid()
ax2.grid()
ax1.set_yscale("log")
ax2.legend()
ax1.set_title("Ring 22")
ax2.set_title("Ring 91")
fig.supxlabel("Distance from inner target (m)")
ax1.set_ylabel("W Density (m-3)")
fig.tight_layout()
fig.show()

# RBS data.
a2path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2path)
a2_shift = 0.01
a2_itf_x = a2["R D (cm)"].values / 100 + a2_shift
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values * 1e15 * 1e4
a2_otf_x = a2["R U (cm)"].values / 100 + a2_shift
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values * 1e15 * 1e4

# Corresponding 3DLIM runs.
root = "/Users/zamperini/Documents/d3d_work/lim_runs/blobtest/167196-a2-tor240-div20"
absfacs = [1.402E13, 5.922E14, 1.144E15, 1.128E15, 1.050E15]
rtip = 2.282
lps = []
for case in tqdm(cases):
    path = root + case + ".nc"
    lps.append(LimPlots.LimPlots(path))

xs_otf = []; ys_otf = []
xs_itf = []; ys_itf = []
absfac_mod = 25 * 5
for i in range(0, len(lps)):
    cent = lps[i].centerline(showplot=False)
    xs_otf.append(cent["x1"] + rtip)
    ys_otf.append(cent["y1"] * absfacs[i] * absfac_mod)
    xs_itf.append(cent["x2"] + rtip)
    ys_itf.append(cent["y2"] * absfacs[i] * absfac_mod)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
ax1.scatter(a2_itf_x, a2_itf_y, color="tab:purple", marker="*", s=100, ec="k", zorder=10)
ax2.scatter(a2_otf_x, a2_otf_y, color="tab:red", marker="*", s=100, ec="k", zorder=10)
for i in range(0, len(lps)):
    ax1.plot(xs_itf[i], ys_itf[i], label=corrs[i], zorder=1)
    ax2.plot(xs_otf[i], ys_otf[i], label=corrs[i], zorder=1)
fig.supxlabel("R (m)")
ax1.set_yscale("log")
ax1.set_ylabel("W Areal Density (m-2)")
ax1.set_title("ITF")
ax2.set_title("OTF")
ax1.set_xlim([2.28, 2.35])
ax1.set_ylim([0, 1e19])
ax2.legend()
fig.tight_layout()
fig.show()
