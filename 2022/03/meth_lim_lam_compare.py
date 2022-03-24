# Script to compare the measured LAMS data to the 3DLIM data.
import LimPlots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
from importlib import reload


# Paths.
shot = 184267
ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-008a.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-001.nc"
#xlpath = "/Users/zamperini/My Drive/Research/Documents/2021/11/test_div_inj.xlsx"
xlpath = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"

# Load in everything.
# 184267: R = ITF, L = OTF
# 184527: R = OTF, L = ITF
reload(LimPlots)
lp = LimPlots.LimPlots(ncpath)
if shot == 184527:
    #itf = pd.read_excel(xlpath, sheet_name="MCP4 L1 Data")
    #otf = pd.read_excel(xlpath, sheet_name="MCP4 R1 Data")
    df = pd.read_excel(xlpath, sheet_name="For Export")
    lams_itf_x = df["ml04_loc"]
    lams_itf_y = df["ml04_excess_c13"]
    lams_otf_x = df["mr04_loc"]
    lams_otf_y = df["mr04_excess_c13"]
elif shot == 184267:
    df = pd.read_excel(xlpath, sheet_name="For Export")
    lams_itf_x = df["mr21_loc"]
    lams_itf_y = df["mr21_excess_c13"]
    lams_otf_x = df["ml21_loc"]
    lams_otf_y = df["ml21_excess_c13"]

# Remove data before zero, smoothing, background subtraction.
mask1 = lams_itf_x > 0
mask2 = lams_otf_x > 0
lams_itf_x = lams_itf_x[mask1]
lams_itf_y = lams_itf_y[mask1]
lams_otf_x = lams_otf_x[mask2]
lams_otf_y = lams_otf_y[mask2]

lams_itf_ys = savgol_filter(lams_itf_y, 51, 2)
lams_otf_ys = savgol_filter(lams_otf_y, 51, 2)

min_itfy = lams_itf_ys.min()
min_otfy = lams_otf_ys.min()
lams_itf_y = lams_itf_y - min_itfy
lams_otf_y = lams_otf_y - min_otfy
lams_itf_ys = lams_itf_ys - min_itfy
lams_otf_ys = lams_otf_ys - min_otfy

# Run our stupid little simulation for a reerosion estimate.
launch_energy      = 10
mass_ion           = 12   # C = 12, Si = 28
charge_ion         = 1
delta_t            = 1e-7
ff_mod             = 1
cdf_opt            = 3
cdf_te_cutoff      = 6
gauss_width        = 500
nparts             = 1000
rad_vel            = None
cdf_exp_falloff    = 0.01
reod1 = lp.sim_reerosion(1, launch_energy=launch_energy, charge_ion=charge_ion,
  ti_mult=1, delta_t=delta_t, ff_mod=ff_mod, rad_vel=rad_vel,
  cdf_opt=cdf_opt, gauss_width=gauss_width, strength=0.0, nparts=nparts,
  cdf_te_cutoff=cdf_te_cutoff, cdf_opt_range=[0.0, 0.1],
  cdf_exp_falloff=cdf_exp_falloff, nruns=1)
reod2 = lp.sim_reerosion(2, launch_energy=launch_energy, charge_ion=charge_ion,
  ti_mult=1, delta_t=delta_t, ff_mod=ff_mod, rad_vel=rad_vel, te_mult=1,
  cdf_opt=cdf_opt, gauss_width=gauss_width, strength=0.0, nparts=nparts,
  cdf_te_cutoff=cdf_te_cutoff, cdf_opt_range=[0.0, 0.1],
  cdf_exp_falloff=cdf_exp_falloff, nruns=1)

# Pull out the data we want to plot.
#lams_itf_x     = itf["Axial Location [mm]"] / 10  # mm to cm
#lams_itf_y     = itf["13C Excess"]
#lams_itf_y_err = itf["13C Excess error"]
#lams_otf_x     = otf["Axial Location (mm)"] / 10
#lams_otf_y     = otf["13C Excess"]
#lams_otf_y_err = otf["13C Excess error"]

# 184267: 1 = ITF, 2 = OTF
# 184527: 1 = OTF, 2 = ITF
lp_data = lp.centerline()
lim_itf_x = lp_data["x2"] * 100  # m to cm
lim_itf_y = lp_data["y2"]
lim_otf_x = lp_data["x1"] * 100
lim_otf_y = lp_data["y1"]

lim_itf_ys = savgol_filter(lim_itf_y, 11, 2)
lim_otf_ys = savgol_filter(lim_otf_y, 11, 2)

show_lim = True

fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(10, 4))
ax2 = ax1.twinx()
ax4 = ax3.twinx()
ax1.plot(lams_itf_x, lams_itf_y, color="tab:red", alpha=0.3)
ax1.plot(lams_itf_x, lams_itf_ys, color="tab:red", label="ITF", lw=3)
if show_lim:
    ax2.plot(lim_itf_x, lim_itf_y, color="k", lw=3)
    ax2.plot(lim_itf_x, lim_itf_y, color="tab:red", lw=2)
ax2.plot(reod2["dep_x"][1:]*100, reod2["net_y"][1:], color="k", lw=2, linestyle="--")
ax3.plot(lams_otf_x, lams_otf_y, color="tab:purple", alpha=0.3)
ax3.plot(lams_otf_x, lams_otf_ys, color="tab:purple", label="OTF", lw=3)
if show_lim:
    ax4.plot(lim_otf_x, lim_otf_y, color="k", lw=3)
    ax4.plot(lim_otf_x, lim_otf_y, color="tab:purple", lw=2)
ax4.plot(reod1["dep_x"][:-1]*100, reod1["net_y"][:-1], color="k", lw=2, linestyle="--")
ax1.set_ylim([0, 1000])
ax3.set_ylim([0, 1000])
ax2.set_ylim([0, 20])
ax4.set_ylim([0, 20])
ax1.legend()
ax3.legend()
ax1.set_xlabel("Distance along probe (cm)")
ax3.set_xlabel("Distance along probe (cm)")
ax1.set_ylabel("LAMS Counts")
ax3.set_ylabel("LAMS Counts")
ax4.set_ylabel("3DLIM Counts")
fig.tight_layout()
fig.show()
