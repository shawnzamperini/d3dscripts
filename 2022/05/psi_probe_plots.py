# Script to make plots of the deposition profiles for PSI.
import pandas as pd
import matplotlib.pyplot as plt
import LimPlots
import numpy as np
from scipy.signal import savgol_filter
from matplotlib.patches import Rectangle


# Inputs
xlpath = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"
ncpath6 = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006.nc"
ncpath6a = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006a.nc"
ncpath6b = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006b.nc"
ncpath6c = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006c.nc"
ncpath6d = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006d.nc"
df = pd.read_excel(xlpath, sheet_name="For Export")
lams_itf_x = df["mr21_loc"]
lams_itf_y = df["mr21_excess_c13"]
lams_otf_x = df["ml21_loc"]
lams_otf_y = df["ml21_excess_c13"]
lim_ylim = [0, 1e16]
lam_ylim = [0, 1200]
absfac = 2.330e+14 * 2.3
show_lim = True
itf_only = True
lim_shift = -0.4

def load_lim_dep(ncpath):
    lp = LimPlots.LimPlots(ncpath)
    lp_data = lp.centerline(showplot=False)
    lim_itf_x = lp_data["x1"] * 100  # m to cm
    lim_itf_y = lp_data["y1"] * absfac
    lim_otf_x = lp_data["x2"] * 100
    lim_otf_y = lp_data["y2"] * absfac

    lim_itf_ys = savgol_filter(lim_itf_y, 11, 2)
    lim_otf_ys = savgol_filter(lim_otf_y, 11, 2)

    # Trim the data near the tip.
    cutx = 0.0
    mask1 = lim_itf_x > cutx
    mask2 = lim_otf_x > cutx
    lim_itf_x = lim_itf_x[mask1]
    lim_itf_y = lim_itf_y[mask1]
    lim_itf_ys = lim_itf_ys[mask1]
    lim_otf_x = lim_otf_x[mask2]
    lim_otf_y = lim_otf_y[mask2]
    lim_otf_ys = lim_otf_ys[mask2]
    return {"itfx": lim_itf_x, "itfy":lim_itf_y, "itfys":lim_itf_ys,
        "otfx":lim_otf_x, "otfy":lim_otf_y, "otfys": lim_otf_ys}

lim6 = load_lim_dep(ncpath6)
lim6a = load_lim_dep(ncpath6a)
lim6b = load_lim_dep(ncpath6b)
lim6c = load_lim_dep(ncpath6c)
lim6d = load_lim_dep(ncpath6d)

# First a plot of the silicon to show net erosion.
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(6,4))

ax1.axhline(3.3e5, color="k", linestyle="--", lw=2)
ax2.axhline(4.42e5, color="k", linestyle="--", lw=2)
ax1.plot(lams_itf_x, savgol_filter(df["mr21_tot_si"], 91, 2), label="ITF", lw=3, color="tab:purple")
ax2.plot(lams_otf_x, savgol_filter(df["ml21_tot_si"], 91, 2), label="OTF", lw=3, color="tab:purple")

ax1.set_xlim([0, 10])
ax1.set_ylim([0, 1.3e6])
ax2.set_xlabel("Distance along probe (cm)", fontsize=16)
fig.supylabel("Si Counts (Arbitrary)", fontsize=16)
fig.tight_layout()
fig.show()

# Remove data before zero, smoothing, background subtraction.
mask1 = lams_itf_x > 0
mask2 = lams_otf_x > 0
lams_itf_x = lams_itf_x[mask1]
lams_itf_y = lams_itf_y[mask1]
lams_otf_x = lams_otf_x[mask2]
lams_otf_y = lams_otf_y[mask2]

# Convert LAMS counts to 1e17 atoms/cm2, and then to atoms/m2.
lams_itf_y = (lams_itf_y + 346.05) / 11942
lams_otf_y = (lams_otf_y + 346.05) / 11942
#lams_itf_y = lams_itf_y * 1e4 * 1e17
#lams_otf_y = lams_otf_y * 1e4 * 1e17
lams_itf_y = lams_itf_y * 1e17
lams_otf_y = lams_otf_y * 1e17

lams_itf_ys = savgol_filter(lams_itf_y, 91, 2)
lams_otf_ys = savgol_filter(lams_otf_y, 91, 2)
min_itfy = lams_itf_ys.min()
min_otfy = lams_otf_ys.min()
lams_itf_y = lams_itf_y - min_itfy
lams_otf_y = lams_otf_y - min_otfy
lams_itf_ys = lams_itf_ys - min_itfy
lams_otf_ys = lams_otf_ys - min_otfy

# Now a plot of the 13C deposition.
if itf_only:
    fig, ax1 = plt.subplots(figsize=(6, 5))
else:
    fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(6, 5), sharex=True, sharey=True)

ax1.plot(lams_itf_x, lams_itf_y, color="k", alpha=0.3)
ax1.plot(lams_itf_x, lams_itf_ys, color="k", label="ITF", lw=3)
if show_lim:

    mask1 = lim6["itfx"] > 0.95

    # Include curves from a convective velocity scan.
    ax1.plot(lim6b["itfx"][mask1]+lim_shift, lim6b["itfys"][mask1], color="k", lw=4)
    ax1.plot(lim6b["itfx"][mask1]+lim_shift, lim6b["itfys"][mask1], color="tab:purple", lw=3)
    ax1.plot(lim6d["itfx"][mask1]+lim_shift, lim6d["itfys"][mask1], color="k", lw=4)
    ax1.plot(lim6d["itfx"][mask1]+lim_shift, lim6d["itfys"][mask1], color="tab:green", lw=3)


    ax1.plot(lim6["itfx"][mask1]+lim_shift, lim6["itfys"][mask1], color="k", lw=4)
    ax1.plot(lim6["itfx"][mask1]+lim_shift, lim6["itfys"][mask1], color="tab:red", lw=3)

if not itf_only:
    ax3.plot(lams_otf_x, lams_otf_y, color="tab:red", alpha=0.3)
    ax3.plot(lams_otf_x, lams_otf_ys, color="tab:red", label="OTF", lw=3)

    # OTF LAMS was re-erosion before ~ 4cm, so cover that region to avoid confusion.
    rect = Rectangle((0, 0), 4, 2e16, color="grey", zorder=10)
    ax3.add_patch(rect)

    if show_lim:
        mask2 = lim6["otfx"] > 0.95
        #ax3.plot(lim6a["otfx"][mask1], lim6a["otfy"][mask1], color="k", lw=2, linestyle="--")
        #ax3.plot(lim6d["otfx"][mask1], lim6d["otfy"][mask1], color="k", lw=2, linestyle="--")
        ax3.plot(lim6["otfx"][mask2], lim6["otfy"][mask2], color="k", lw=3, zorder=20)
        ax3.plot(lim6["otfx"][mask2], lim6["otfy"][mask2], color="tab:red", lw=2, zorder=20)

    ax3.tick_params(labelsize=14)
    ax3.legend(fontsize=14)
    ax3.set_xlabel("Distance along probe (cm)", fontsize=16)
    fig.supylabel(r"$\mathdefault{^{13}C\ Areal\ Density\ (atoms/cm^2)}$", fontsize=16)

else:
    ax1.set_xlabel("Distance along probe (cm)", fontsize=16)
    ax1.set_ylabel(r"$\mathdefault{^{13}C\ Areal\ Density\ (atoms/cm^2)}$", fontsize=16)

ax1.tick_params(labelsize=14)
ax1.legend(fontsize=14)
ax1.set_xlim([0, 6])
ax1.set_ylim(lim_ylim)
fig.tight_layout()
fig.show()
