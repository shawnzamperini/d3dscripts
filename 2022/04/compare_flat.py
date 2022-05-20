# This script compares a full blown coupled DIVIMP-3DLIM run with a comparable
# 3DLIM run with just a flat impurity distribution. See
# ../03/meth_lim_lam_compare.py for a more complicated script that this one
# is born out of.
import LimPlots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from importlib import reload
from scipy.signal import savgol_filter
from matplotlib.lines import Line2D


# Inputs
xlpath = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"
xlim = [0, 6]
ylim = [0, 2e16]
ncpath527 = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-020.nc"
ncpath267 = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006.nc"
ncpath527_flat = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-021.nc"
ncpath267_flat = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-007.nc"


reload(LimPlots)

# 184267: 1 = ITF, 2 = OTF
# 184527: 1 = OTF, 2 = ITF
def load_data(ncpath, shot):

    lp = LimPlots.LimPlots(ncpath)
    lp_data = lp.centerline(showplot=False)

    if shot == 184527:

        df = pd.read_excel(xlpath, sheet_name="For Export")
        lams_itf_x = df["ml04_loc"]
        lams_itf_y = df["ml04_excess_c13"]
        lams_otf_x = df["mr04_loc"]
        lams_otf_y = df["mr04_excess_c13"]

        # ABSFAC from 015 DIVIMP simulation. Need to multiply by 2.3 since UOB
        # was open for 2.3 seconds.
        absfac = 4.749e14 * 2.3
        lim_itf_x = lp_data["x2"] * 100  # m to cm
        lim_itf_y = lp_data["y2"] * absfac
        lim_otf_x = lp_data["x1"] * 100
        lim_otf_y = lp_data["y1"] * absfac

        # Load NRA data. Only available for this shot.
        df = pd.read_excel(xlpath, sheet_name="NRA")
        nra_itf_c12 = df["ml04_c12 (1e17 cm-2)"] * 1e17
        nra_itf_c13 = df["ml04_c13 (1e17 cm-2)"] * 1e17
        nra_itf_x = df["ml04_x (mm)"] / 10
        nra_itf_y = nra_itf_c13 - (0.01/0.99) * nra_itf_c12
        nra_otf_c12 = df["mr04_c12 (1e17 cm-2)"] * 1e17
        nra_otf_c13 = df["mr04_c13 (1e17 cm-2)"] * 1e17
        nra_otf_x = df["mr04_x (mm)"] / 10
        nra_otf_y = nra_otf_c13 - (0.01/0.99) * nra_otf_c12

    elif shot == 184267:

        df = pd.read_excel(xlpath, sheet_name="For Export")
        lams_itf_x = df["mr21_loc"]
        lams_itf_y = df["mr21_excess_c13"]
        lams_otf_x = df["ml21_loc"]
        lams_otf_y = df["ml21_excess_c13"]

        # ABSFAC from 015 DIVIMP simulation. Need to multiply by 2.3 since UOB
        # was open for 2.3 seconds.
        absfac = 2.330e+14 * 2.3
        lim_itf_x = lp_data["x1"] * 100  # m to cm
        lim_itf_y = lp_data["y1"] * absfac
        lim_otf_x = lp_data["x2"] * 100
        lim_otf_y = lp_data["y2"] * absfac

    # Remove data before zero, smoothing, background subtraction.
    mask1 = lams_itf_x > 0
    mask2 = lams_otf_x > 0
    lams_itf_x = lams_itf_x[mask1]
    lams_itf_y = lams_itf_y[mask1]
    lams_otf_x = lams_otf_x[mask2]
    lams_otf_y = lams_otf_y[mask2]

    # Convert LAMS counts to 1e17 atoms/cm2.
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

    lim_itf_ys = savgol_filter(lim_itf_y, 11, 2)
    lim_otf_ys = savgol_filter(lim_otf_y, 11, 2)

    if shot == 184267:
        return {"lim_itf_x":lim_itf_x, "lim_itf_y":lim_itf_y,
            "lim_otf_x":lim_otf_x, "lim_otf_y":lim_otf_y,
            "lams_itf_x":lams_itf_x, "lams_itf_ys":lams_itf_ys,
            "lams_otf_x":lams_otf_x, "lams_otf_ys":lams_otf_ys}
    elif shot == 184527:
        return {"lim_itf_x":lim_itf_x, "lim_itf_y":lim_itf_y,
            "lim_otf_x":lim_otf_x, "lim_otf_y":lim_otf_y,
            "lams_itf_x":lams_itf_x, "lams_itf_ys":lams_itf_ys,
            "lams_otf_x":lams_otf_x, "lams_otf_ys":lams_otf_ys,
            "nra_itf_x":nra_itf_x, "nra_itf_y":nra_itf_y,
            "nra_otf_x":nra_otf_x, "nra_otf_y":nra_otf_y}

# Load in all the data.
d527 = load_data(ncpath527, 184527)
d267 = load_data(ncpath267, 184267)
d527f = load_data(ncpath527_flat, 184527)
d267f = load_data(ncpath267_flat, 184267)


lw = 3
plt.rcParams["font.family"] = "Century Gothic"
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 5), sharex=True, sharey=True)

# One large Axes to hold labels.
ax = fig.add_subplot(111, frameon=False)
ax.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
ax.set_xlabel("Distance along probe (cm)", fontsize=16)
ax.set_ylabel(r"$\mathdefault{^{13}C\ Areal\ Density\ (atoms/cm^2)}$", fontsize=16)

# 184527 plots on top row.
ax1.plot(d527["lams_itf_x"], d527["lams_itf_ys"], color="tab:purple", lw=lw, label="LAMS")
ax1.plot(d527["lim_itf_x"], d527["lim_itf_y"], color="k", lw=lw+1)
ax1.plot(d527["lim_itf_x"], d527["lim_itf_y"], color="tab:purple", lw=lw, label="3DLIM")
ax1.plot(d527f["lim_itf_x"], d527f["lim_itf_y"], color="k", lw=lw-1)

ax2.plot(d527["lams_otf_x"], d527["lams_otf_ys"], color="tab:purple", lw=lw, label="LAMS")
ax2.plot(d527["lim_otf_x"], d527["lim_otf_y"], color="k", lw=lw+1)
ax2.plot(d527["lim_otf_x"], d527["lim_otf_y"], color="tab:purple", lw=lw, label="3DLIM")
ax2.plot(d527f["lim_otf_x"], d527f["lim_otf_y"], color="k", lw=lw-1)

ax3.plot(d267["lams_itf_x"], d267["lams_itf_ys"], color="tab:red", lw=lw, label="LAMS")
ax3.plot(d267["lim_itf_x"], d267["lim_itf_y"], color="k", lw=lw+1)
ax3.plot(d267["lim_itf_x"], d267["lim_itf_y"], color="tab:red", lw=lw, label="3DLIM")
ax3.plot(d267f["lim_itf_x"], d267f["lim_itf_y"], color="k", lw=lw-1)

ax4.plot(d267["lams_otf_x"], d267["lams_otf_ys"], color="tab:red", lw=lw, label="LAMS")
ax4.plot(d267["lim_otf_x"], d267["lim_otf_y"], color="k", lw=lw+1)
ax4.plot(d267["lim_otf_x"], d267["lim_otf_y"], color="tab:red", lw=lw, label="3DLIM")
ax4.plot(d267f["lim_otf_x"], d267f["lim_otf_y"], color="k", lw=lw-1)

ax.set_title("ITF" + " "*50 + "OTF", fontsize=16)
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)

fig.tight_layout()
fig.show()
