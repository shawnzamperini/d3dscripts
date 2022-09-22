# Quick script to just plot some example RCP data.
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from scipy.optimize import curve_fit


rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
rcp = pd.read_excel(rcp_path, sheet_name="MP184267_1")

x = rcp["Dist. from sep."].values * 100
y = rcp["ne (1e18 m-3)"].values * 1e18
qp = rcp["q_par (MW/m2)"].values

plt.rcParams['font.family'] = 'Century Gothic'
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(x, y, color="k", lw=3)
ax.set_xlabel("Distance from separatrix (cm)", fontsize=14)
ax.set_ylabel(r"$\mathdefault{n_e\ (m^{-3})}$", fontsize=14)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_yscale("log")
ax.set_ylim(1e17, 1e19)
ax.tick_params(labelsize=12, which="both")
#ax.grid(axis="y", which="both")
#ax.grid(axis="x", which="major")
fig.tight_layout()
fig.show()

# Exponential fit.
def exp_fit(x, a, b):
    return a * np.exp(-b*x)
mask = np.logical_and(x > 6, x<9)
popt, pcov = curve_fit(exp_fit, x[mask], qp[mask])

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(x[mask], exp_fit(x[mask], *popt), color="k", lw=3, zorder=1)
ax.scatter(x, qp, color="tab:red", s=45, marker="^", edgecolors="k", linewidths=1, zorder=5)
ax.set_xlabel("Distance from separatrix (cm)", fontsize=14)
ax.set_ylabel(r"$\mathdefault{q_{\|\|}\ (MW/m^{2})}$", fontsize=14)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_yscale("log")
ax.set_xlim([2, 9])
#ax.set_ylim(1e17, 1e19)
ax.tick_params(labelsize=12, which="both")
#ax.grid(axis="y", which="both")
#ax.grid(axis="x", which="major")
fig.tight_layout()
fig.show()

# Estimate the SOL collionality parameter.
lconn_path = "/Users/zamperini/My Drive/Research/Documents/2022/03/lconns_184527_mimes.xlsx"
df = pd.read_excel(lconn_path)
lconn_r = df["R (m)"].values * 100
lconn = df["Lconn (km)"].values * 1000
r = rcp["R (cm)"].values
ne = y
te = rcp["Te (eV)"].values
ti_near = te[-1]
ti_far = te[0] * 4

f_te = interp1d(r, te)
f_ne = interp1d(r, ne)
f_lconn = interp1d(lconn_r, lconn)
f_ti = interp1d([r[-1], r[0]], [ti_near, ti_far])

nu_sol = 1e-16 * f_ne(r) * f_lconn(r) / np.square(f_ti(r))

rsep = 226
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(r-rsep, nu_sol, color="k", lw=3)
ax.set_xlabel("R-Rsep (cm)", fontsize=16)
ax.set_ylabel(r"$\mathdefault{\nu_{SOL}^*}$", fontsize=16)
ax.tick_params(labelsize=14)
ax.set_yscale("log")
ax.grid(which="both", alpha=0.4)
fig.tight_layout()
fig.show()
