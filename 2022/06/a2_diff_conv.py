# Script to compare a diffusive simulation from DIVIMP to a convective one
# using the DIVIMP CP model.
import oedge_plots
import LimPlots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


exp_time = 4

# Diffusive model.
absfac_mod = 1.0
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-diffusion.nc"
op = oedge_plots.OedgePlots(ncpath)
diff_cp = op.collector_probe(2.27, -0.188, 2.37, -0.188, showplot=False, numlocs=30)
diff_cp_itf_x = np.array(diff_cp["r"])
diff_cp_itf_y = np.array(diff_cp["flux1"]) * exp_time * absfac_mod
diff_cp_otf_x = np.array(diff_cp["r"])
diff_cp_otf_y = np.array(diff_cp["flux2"]) * exp_time * absfac_mod
mask = ~np.isnan(diff_cp_itf_y)
diff_cp_itf_x = diff_cp_itf_x[mask]
diff_cp_itf_y = diff_cp_itf_y[mask]

# Convective model.
absfac_mod = 0.25
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div11.nc"
op = oedge_plots.OedgePlots(ncpath)
conv_cp = op.collector_probe(2.27, -0.188, 2.37, -0.188, showplot=False, numlocs=30)
conv_cp_itf_x = np.array(conv_cp["r"])
conv_cp_itf_y = np.array(conv_cp["flux1"]) * exp_time * absfac_mod
conv_cp_otf_x = np.array(conv_cp["r"])
conv_cp_otf_y = np.array(conv_cp["flux2"]) * exp_time * absfac_mod
mask = ~np.isnan(conv_cp_itf_y)
conv_cp_itf_x = conv_cp_itf_x[mask]
conv_cp_itf_y = conv_cp_itf_y[mask]

# Load 3DLIM run. Absfac * shots * time of shot * fudge factor
absfac = 1.647e14 * 25 * 5 * 10
rtip = 2.282
ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/blobtest/167196-a2-tor240_45-blob.nc"
lp = LimPlots.LimPlots(ncpath)
cent = lp.centerline(showplot=False)
lim_otfx = cent["x1"] + rtip
lim_otfy = cent["y1"] * absfac
lim_itfx = cent["x2"] + rtip
lim_itfy = cent["y2"] * absfac

# Get A2 data from sheet.
a2path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2path)
a2_shift = 0.01
a2_itf_x = a2["R D (cm)"].values / 100 + a2_shift
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values * 1e15 * 1e4
a2_otf_x = a2["R U (cm)"].values / 100 + a2_shift
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values * 1e15 * 1e4

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(8, 4))

lw = 3
ax1.plot(diff_cp_itf_x, diff_cp_itf_y, color="tab:purple", linestyle="--", lw=lw)
ax2.plot(diff_cp_otf_x, diff_cp_otf_y, color="tab:red", linestyle="--", lw=lw, label="D=0.3 m2/s")
ax1.plot(conv_cp_itf_x, conv_cp_itf_y, color="tab:purple", lw=lw)
ax2.plot(conv_cp_otf_x, conv_cp_otf_y, color="tab:red", lw=lw, label="vmean=500 m/s")

#ax11 = ax1.twinx()
#ax22 = ax2.twinx()
ax1.plot(lim_itfx, lim_itfy, color="k")
ax2.plot(lim_otfx, lim_otfy, color="k", label="3DLIM")

ax1.scatter(a2_itf_x, a2_itf_y, color="tab:purple", marker="*", s=75, ec="k")
ax2.scatter(a2_otf_x, a2_otf_y, color="tab:red", marker="*", s=75, ec="k")
fig.supxlabel("R (m)")
ax1.set_ylabel("W Areal Density (m-2)")
ax1.set_title("ITF")
ax2.set_title("OTF")
ax1.set_xlim([2.28, 2.35])
ax1.set_ylim([0, 1e19])
ax2.legend()
fig.tight_layout()
fig.show()
