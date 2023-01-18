import pandas as pd
from scipy.interpolate import interp1d, splrep, BSpline, griddata
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt, savgol_filter


# Fudge factors.
mp_psin_shift = -0.04
ts_ne_mult = 1.12
mp_ne_mult = 0.75

# Load data and gfile.
xlpath = "/Users/zamperini/Documents/d3d_work/files/solps_data_eric/167196_data.xlsx"
xl_lp = pd.read_excel(xlpath, sheet_name="Outer Target")
xl_ts = pd.read_excel(xlpath, sheet_name="Thomson Scattering")
xl_mp = pd.read_excel(xlpath, sheet_name="RCP Data")
xl_bl = pd.read_excel(xlpath, sheet_name="Blob Data")
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)


# Spline fit.
def spline_fit(x, y_in, knot_numbers=8, s=1):
    maxy = y_in.max()
    y = y_in / maxy
    x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
    q_knots = np.quantile(x, x_new)
    t,c,k = splrep(x, y, t=q_knots, s=s)
    yfit = BSpline(t, c, k)(x)
    return yfit * maxy

# Just median filter for LP data.
lp_psin = xl_lp["Psin"].values
lp_ne = xl_lp["ne (m-3)"].values
lp_te = xl_lp["Te (eV)"].values
lp_qpar = xl_lp["qpar (W/m2)"].values
#lp_ne_fit = spline_fit(lp_psin, lp_ne)
lp_ne_fit = medfilt(lp_ne, 25)
#lp_te_fit = spline_fit(lp_psin, lp_te)
lp_te_fit = medfilt(lp_te, 15)
lp_qpar_fit = medfilt(lp_qpar, 21)
#lp_qpar_fit = savgol_filter(lp_qpar, 15, 2)

# TS data just average the fits together.
ts_mean = xl_ts.groupby("Psin").mean()
ts_psin = ts_mean.index
ts_ne = ts_mean["ne (m-3)"].values * ts_ne_mult
ts_te = ts_mean["Te (eV)"].values

# RCP do a spline.
mp_psin = xl_mp["Psin"].values
sort_idx = np.argsort(mp_psin)
mp_psin = mp_psin[sort_idx]
mp_ne = xl_mp["Ne(E18 m-3)"].values[sort_idx] * 1e18 * mp_ne_mult
mp_te = xl_mp["Te(eV)"].values[sort_idx]
mp_r = xl_mp["R(cm)"].values[sort_idx] / 100
mp_ne_fit = spline_fit(mp_psin, mp_ne)
mp_te_fit = spline_fit(mp_psin, mp_te)
mp_psin = mp_psin + mp_psin_shift

# For the final midplane profiles we will do a hybrid of the TS and RCP data. The SOL will be only the probe data,
# while the core will use TS. In the small gap between the two the data is linearly interpolated.
mid_psin = np.linspace(0, 1.35, 250)
#f_ne = interp1d(ts_psin, ts_ne, fill_value="extrapolate")
tmp_te = ts_te[ts_psin <= 1.0]
tmp_ne = ts_ne[ts_psin <= 1.0]
tmp_psin = ts_psin[ts_psin <= 1.0]
tmp_te = np.append(tmp_te, mp_te_fit[mp_psin > 1.0])
tmp_ne = np.append(tmp_ne, mp_ne_fit[mp_psin > 1.0])
tmp_psin = np.append(tmp_psin, mp_psin[mp_psin > 1.0])
f_te = interp1d(tmp_psin, tmp_te, fill_value="extrapolate")
f_ne = interp1d(tmp_psin, tmp_ne, fill_value="extrapolate")

# Use the gfile to interpolate the psin values to an R-Rsep @ OMP coordinate.
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
R_axis = gfile['RMAXIS']
Z_axis = gfile["ZMAXIS"]
Rs_trunc = Rs > R_axis
print("Making interpolation function...")
mid_romp = griddata((psin[Rs_trunc].flatten(), Zs[Rs_trunc].flatten()), Rs[Rs_trunc].flatten(),
                    list(zip(mid_psin, np.full(len(mid_psin), Z_axis))))
rsepomp = griddata((psin[Rs_trunc].flatten(), Zs[Rs_trunc].flatten()), Rs[Rs_trunc].flatten(), (1.0, Z_axis))
mid_rmrsomp = mid_romp - rsepomp

# Some extra blob columns to add.
f_psin = interp1d(mp_r, mp_psin)
blob_r = xl_bl["R(cm)"].values / 100
blob_psin = f_psin(blob_r)
blob_romp = griddata((psin[Rs_trunc].flatten(), Zs[Rs_trunc].flatten()), Rs[Rs_trunc].flatten(),
                     list(zip(blob_psin, np.full(len(blob_psin), Z_axis))))
blob_rmrsomp = blob_romp - rsepomp

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 7))

ax1.axvline(1.0, color="k")
ax2.axvline(1.0, color="k")
ax3.axvline(1.0, color="k")
ax4.axvline(1.0, color="k")

ax1.scatter(lp_psin, lp_ne, marker="^", edgecolors="k", color="tab:red", s=15, alpha=0.6)
ax1.plot(lp_psin, lp_ne_fit, color="k", lw=3)
ax1.plot(lp_psin, lp_ne_fit, color="tab:red", lw=2)
ax1.set_ylabel("ne (m-3)")

ax3.scatter(lp_psin, lp_te, marker="^", edgecolors="k", color="tab:red", s=15, alpha=0.6)
ax3.plot(lp_psin, lp_te_fit, color="k", lw=3)
ax3.plot(lp_psin, lp_te_fit, color="tab:red", lw=2)
ax3.set_ylabel("Te (eV)")

ax2.scatter(ts_psin, ts_ne, marker="^", edgecolors="k", color="tab:red", s=15, alpha=0.6)
ax2.set_xlabel("Psin")
ax2.set_ylim([0, 2e19])
ax2.scatter(mp_psin, mp_ne, marker="^", edgecolors="k", color="tab:purple", s=15, alpha=0.6)
ax2.plot(mp_psin, mp_ne_fit, color="k", lw=3)
ax2.plot(mp_psin, mp_ne_fit, color="tab:purple", lw=2)
ax2.plot(mid_psin, f_ne(mid_psin), color="k")
ax2.set_xlim([0.98, None])

ax4.scatter(ts_psin, ts_te, marker="^", edgecolors="k", color="tab:red", s=15, alpha=0.6)
ax4.set_xlabel("Psin")
ax4.set_ylim([0, 100])
ax4.scatter(mp_psin, mp_te, marker="^", edgecolors="k", color="tab:purple", s=15, alpha=0.6)
ax4.plot(mp_psin, mp_te_fit, color="k", lw=3)
ax4.plot(mp_psin, mp_te_fit, color="tab:purple", lw=2)
ax4.plot(mid_psin, f_te(mid_psin), color="k")
ax4.set_xlim([0.98, None])

fig.tight_layout()
fig.show()

# One extra plot of the heatflux.
fig, ax = plt.subplots(figsize=(5, 4))
ax.axvline(1.0, color="k")
ax.scatter(lp_psin, lp_qpar, marker="^", edgecolors="k", color="tab:red", s=15, alpha=0.6)
ax.plot(lp_psin, lp_qpar_fit, color="k", lw=3)
ax.plot(lp_psin, lp_qpar_fit, color="tab:red", lw=2)
ax.set_ylabel("qpar (W/m2)")
fig.tight_layout()
fig.show()