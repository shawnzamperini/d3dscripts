import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import LimPlots
from importlib import reload


shot = 184527
ti_mult = 1
exp_time = 5

rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
lams_path = "/Users/zamperini/My Drive/Research/Data/lams_data/methane_lams_master.xlsx"
if shot == 184527:
    rcp1 = pd.read_excel(rcp_path, sheet_name="MP184527_1")
    rcp2 = pd.read_excel(rcp_path, sheet_name="MP184527_2")
    rcp_r  = rcp2["R (cm)"].values[:-7]
    te = rcp2["Te (eV)"].values[:-7]
    ne = rcp2["ne (1e18 m-3)"].values[:-7] * 1e18
    te[rcp_r > 235] = 3
    ne[rcp_r > 235] = 3e17
    ptip = 2.318

    df = pd.read_excel(lams_path, sheet_name="For Export")
    lams_itf_x = df["ml04_loc"].values
    lams_itf_y = df["ml04_excess_c13"].values
    lams_otf_x = df["mr04_loc"].values
    lams_otf_y = df["mr04_excess_c13"].values

    # A 3DLIM run.
    ncpath = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-008.nc"
    absfac = 2.102e15  # C13 atoms/m3, from DIVIMP?
    reload(LimPlots)
    lp = LimPlots.LimPlots(ncpath)
    lp_data = lp.centerline()
    lim_itf_x = lp_data["x2"] * 100  # m to cm
    lim_itf_y = lp_data["y2"]
    lim_otf_x = lp_data["x1"] * 100
    lim_otf_y = lp_data["y1"]

# Load in chemical D-->C yield data.
yield_path = "/Users/zamperini/My Drive/Research/Documents/2022/02/mixed_material.xlsx"
yields = pd.read_excel(yield_path, sheet_name="Yields", skiprows=5)
y_e = yields["E.3"].values
y_ch = yields["Ysurf"].values

# Convert LAMS counts to 1e17 atoms/cm2, and then to atoms/m2.
lams_itf_y = (lams_itf_y + 346.05) / 11942
lams_otf_y = (lams_otf_y + 346.05) / 11942
lams_itf_y = lams_itf_y * 1e4 * 1e17
lams_otf_y = lams_otf_y * 1e4 * 1e17

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

# Convert to CP coordinates, smoothing, interpolations.
rcp_cp_x = rcp_r - ptip * 100
te_s = savgol_filter(te, 7, 2)
ne_s = savgol_filter(ne, 7, 2)
f_te = interp1d(rcp_cp_x, te_s, fill_value="extrapolate")
f_ne = interp1d(rcp_cp_x, ne_s, fill_value="extrapolate")
f_y = interp1d(y_e, y_ch)
f_itfy = interp1d(lams_itf_x, lams_itf_ys)
f_otfy = interp1d(lams_otf_x, lams_otf_ys)

# Calculate yields along probe surface.
comx = np.linspace(0, 8, 100)
eimp = 3 * f_te(comx) + 2 * ti_mult * f_te(comx)
ys = f_y(eimp)
mb = 2.0 * 931.49 * 10**6 / ((3*10**8)**2)   # eV * s2 / m2
cs = np.sqrt((f_te(comx) + ti_mult * f_te(comx)) / mb)

# Eroded flux is cs * yield.
flux_ero = f_ne(comx) * cs * ys

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

ax1.scatter(rcp_cp_x, te, s=10, color="r")
ax1.plot(comx, f_te(comx), color="r")
ax1.set_xlabel("Distance along probe (cm)")
ax1.set_ylabel("Te (eV)")

ax2.scatter(rcp_cp_x, ne, s=10, color="g")
ax2.plot(comx, f_ne(comx), color="g")
ax2.set_xlabel("Distance along probe (cm)")
ax2.set_ylabel("ne (m-3)")

ax3.plot(lams_itf_x, lams_itf_y/exp_time, alpha=0.2)
ax3.plot(lams_itf_x, lams_itf_ys/exp_time)
ax3.set_xlabel("Distance along probe (cm)")
ax3.set_ylabel("Net C13 Flux (atoms/m2/s)")

#ax33 = ax3.twinx()
ax3.plot(comx, flux_ero, color="k", linestyle="--")
#ax33.set_ylabel("Eroded Flux")

ax4.plot(lams_otf_x, lams_otf_y/exp_time, alpha=0.2)
ax4.plot(lams_otf_x, lams_otf_ys/exp_time)
ax4.set_xlabel("Distance along probe (cm)")
ax4.set_ylabel("C13 atoms/m2")

#ax44 = ax4.twinx()
ax4.plot(comx, flux_ero, color="k", linestyle="--")
#ax44.set_ylabel("Eroded Flux")

fig.tight_layout()
fig.show()
