# Script to compare the measured LAMS data to the 3DLIM data.
import LimPlots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
from importlib import reload


# Inputs
shot = 184267
show_reerosion = True
lim_sameaxis = True
xlpath = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"

# Load in everything.
# 184267: R = ITF, L = OTF
# 184527: R = OTF, L = ITF
print(shot)
reload(LimPlots)

if shot == 184527:
    ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-020.nc"
    #lim_ylim = [0, 22]
    #lam_ylim = [0, 2000]
    lim_ylim = [0, 1.5e16]
    lam_ylim = [0, None]
    df = pd.read_excel(xlpath, sheet_name="For Export")
    lams_itf_x = df["ml04_loc"]
    lams_itf_y = df["ml04_excess_c13"]
    lams_otf_x = df["mr04_loc"]
    lams_otf_y = df["mr04_excess_c13"]
    linecolor = "tab:purple"

    # ABSFAC from 015 DIVIMP simulation. Need to multiply by 2.3 since UOB was
    # open for 2.3 seconds.
    absfac = 4.749e14 * 2.3

    # Load NRA data.
    df = pd.read_excel(xlpath, sheet_name="NRA")
    nra_itf_c12 = df["ml04_c12 (1e17 cm-2)"] * 1e17
    nra_itf_c13 = df["ml04_c13 (1e17 cm-2)"] * 1e17
    nra_itf_x = df["ml04_x (mm)"] / 10
    nra_itf_y = nra_itf_c13 - (0.01/0.99) * nra_itf_c12
    nra_otf_c12 = df["mr04_c12 (1e17 cm-2)"] * 1e17
    nra_otf_c13 = df["mr04_c13 (1e17 cm-2)"] * 1e17
    nra_otf_x = df["mr04_x (mm)"] / 10
    nra_otf_y = nra_otf_c13 - (0.01/0.99) * nra_otf_c12

    strength1          = 0.0
    strength2          = 0.0
    cdf_opt_range1     = [0.0, 0.1]
    cdf_opt_range2     = [0.0, 0.1]
    launch_energy      = 10
    mass_ion           = 12   # C = 12, Si = 28
    charge_ion         = 1
    delta_t            = 1e-7
    ff_mod             = 1
    cdf_opt            = 3
    cdf_te_cutoff      = 6
    gauss_width        = 100
    nparts             = 1000
    rad_vel            = None
    cdf_exp_falloff    = 0.01

elif shot == 184267:
    ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006.nc"
    ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-008.nc"
    lim_ylim = [0, 1.5e16]
    lam_ylim = [0, 1200]
    df = pd.read_excel(xlpath, sheet_name="For Export")
    lams_itf_x = df["mr21_loc"]
    lams_itf_y = df["mr21_excess_c13"]
    lams_otf_x = df["ml21_loc"]
    lams_otf_y = df["ml21_excess_c13"]
    linecolor = "tab:red"

    # ABSFAC from 015 DIVIMP simulation. Need to multiply by 2.3 since UOB was
    # open for 2.3 seconds.
    #absfac = 2.330e+14 * 2.3
    absfac = 1.0 * 2.3

    strength1          = 0.0  # ITF
    strength2          = 0.0  # OTF
    cdf_opt_range1     = [0.0, 0.04]
    cdf_opt_range2     = [0.0, 0.03]
    launch_energy      = 10
    mass_ion           = 12   # C = 12, Si = 28
    charge_ion         = 1
    delta_t            = 1e-7
    ff_mod             = 1
    cdf_opt            = 2
    cdf_te_cutoff      = 7
    gauss_width        = 100
    nparts             = 10000
    rad_vel            = None
    cdf_exp_falloff    = 0.03

lp = LimPlots.LimPlots(ncpath)

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

# Run our stupid little simulation for a reerosion estimate.
reod1 = lp.sim_reerosion(1, launch_energy=launch_energy, charge_ion=charge_ion,
  ti_mult=1, delta_t=delta_t, ff_mod=ff_mod, rad_vel=rad_vel,
  cdf_opt=cdf_opt, gauss_width=gauss_width, strength=strength1, nparts=nparts,
  cdf_te_cutoff=cdf_te_cutoff, cdf_opt_range=cdf_opt_range1,
  cdf_exp_falloff=cdf_exp_falloff, nruns=1)
reod2 = lp.sim_reerosion(2, launch_energy=launch_energy, charge_ion=charge_ion,
  ti_mult=1, delta_t=delta_t, ff_mod=ff_mod, rad_vel=rad_vel, te_mult=1,
  cdf_opt=cdf_opt, gauss_width=gauss_width, strength=strength2, nparts=nparts,
  cdf_te_cutoff=cdf_te_cutoff, cdf_opt_range=cdf_opt_range2,
  cdf_exp_falloff=cdf_exp_falloff, nruns=1)

# 184267: 1 = ITF, 2 = OTF
# 184527: 1 = OTF, 2 = ITF
lp_data = lp.centerline(showplot=False)
if shot == 184527:
    lim_itf_x = lp_data["x2"] * 100  # m to cm
    lim_itf_y = lp_data["y2"] * absfac
    lim_otf_x = lp_data["x1"] * 100
    lim_otf_y = lp_data["y1"] * absfac
elif shot == 184267:
    lim_itf_x = lp_data["x1"] * 100  # m to cm
    lim_itf_y = lp_data["y1"] * absfac
    lim_otf_x = lp_data["x2"] * 100
    lim_otf_y = lp_data["y2"] * absfac

lim_itf_ys = savgol_filter(lim_itf_y, 11, 2)
lim_otf_ys = savgol_filter(lim_otf_y, 11, 2)

show_lim = True

fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(10, 4))
if lim_sameaxis:
    ax1.plot(lams_itf_x, lams_itf_y, color=linecolor, alpha=0.3)
    ax1.plot(lams_itf_x, lams_itf_ys, color=linecolor, label="ITF", lw=3)
    ax1.plot(lim_itf_x, lim_itf_y, color="k", lw=3)
    ax1.plot(lim_itf_x, lim_itf_y, color=linecolor, lw=2)
    ax3.plot(lams_otf_x, lams_otf_y, color=linecolor, alpha=0.3)
    ax3.plot(lams_otf_x, lams_otf_ys, color=linecolor, label="OTF", lw=3)
    ax3.plot(lim_otf_x, lim_otf_y, color="k", lw=3)
    ax3.plot(lim_otf_x, lim_otf_y, color=linecolor, lw=2)

    # NRA only available for this probe.
    #if shot == 184527:
    #    ax1.scatter(nra_itf_x, nra_itf_y, color="tab:red", edgecolor="k")
    #    ax3.scatter(nra_otf_x, nra_otf_y, color="tab:purple", edgecolor="k")

    ax1.tick_params(labelsize=14)
    ax3.tick_params(labelsize=14)
    ax1.legend(fontsize=14)
    ax3.legend(fontsize=14)
    ax1.set_xlabel("Distance along probe (cm)", fontsize=16)
    ax3.set_xlabel("Distance along probe (cm)", fontsize=16)
    ax1.set_ylabel(r"$\mathdefault{^{13}C\ Areal\ Density\ (atoms/cm^2)}$", fontsize=16)
    #ax3.set_ylabel("C13 Areal Density (atoms/cm2)")
    ax1.set_ylim(lim_ylim)
    ax3.set_ylim(lim_ylim)
else:
    ax2 = ax1.twinx()
    ax4 = ax3.twinx()
    ax1.plot(lams_itf_x, lams_itf_y, color="tab:red", alpha=0.3)
    ax1.plot(lams_itf_x, lams_itf_ys, color="tab:red", label="ITF", lw=3)
    if show_lim:
        ax2.plot(lim_itf_x, lim_itf_y, color="k", lw=3)
        ax2.plot(lim_itf_x, lim_itf_y, color="tab:red", lw=2)

    ax3.plot(lams_otf_x, lams_otf_y, color="tab:purple", alpha=0.3)
    ax3.plot(lams_otf_x, lams_otf_ys, color="tab:purple", label="OTF", lw=3)
    if show_lim:
        ax4.plot(lim_otf_x, lim_otf_y, color="k", lw=3)
        ax4.plot(lim_otf_x, lim_otf_y, color="tab:purple", lw=2)

    if show_reerosion:
        if shot == 184527:
            ax2.plot(reod2["dep_x"][1:]*100, reod2["net_y"][1:], color="k", lw=2, linestyle="--")
            ax4.plot(reod1["dep_x"][:-1]*100, reod1["net_y"][:-1], color="k", lw=2, linestyle="--")
        elif shot == 184267:
            ax4.plot(reod2["dep_x"][1:]*100, reod2["net_y"][1:], color="k", lw=2, linestyle="--")
            ax2.plot(reod1["dep_x"][:-1]*100, reod1["net_y"][:-1], color="k", lw=2, linestyle="--")
    ax1.set_ylim(lam_ylim)
    ax3.set_ylim(lam_ylim)
    ax2.set_ylim(lim_ylim)
    ax4.set_ylim(lim_ylim)
    ax1.set_xlim([0, 10])
    ax3.set_xlim([0, 10])
    ax1.legend()
    ax3.legend()
    ax1.set_xlabel("Distance along probe (cm)")
    ax3.set_xlabel("Distance along probe (cm)")
    ax1.set_ylabel("LAMS Counts")
    ax3.set_ylabel("LAMS Counts")
    ax4.set_ylabel("3DLIM Counts")
fig.tight_layout()
fig.show()
