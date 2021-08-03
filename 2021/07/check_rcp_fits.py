# Helper script to check is an RCP plunge is easily fit with a line or exponential.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit


plunge = "MP187103_1"
df = pd.read_excel("/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx", sheet_name=plunge)

r = df["R (cm)"] / 100
ne = df["ne (1e18 m-3)"]
te = df["Te (eV)"]
vp = (df["Vf1 (V)"] + df["Vf2 (V)"]) / 2 + 3 * te

# Exponential fits.
def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c
exp_popt_ne, exp_pcov_ne = curve_fit(exp_fit, (r - 2.32), ne, p0=(10, 10, 1), maxfev=5000)
exp_popt_te, exp_pcov_te = curve_fit(exp_fit, (r - 2.32), te, p0=(10, 10, 1), maxfev=5000)
exp_popt_vp, exp_pcov_vp = curve_fit(exp_fit, (r - 2.32), vp, p0=(10, 10, 1), maxfev=5000)

# Linear fits.
def lin_fit(x, m, b):
    return m * x + b
lin_popt_ne, lin_pcov_ne = curve_fit(lin_fit, r, ne)
lin_popt_te, lin_pcov_te = curve_fit(lin_fit, r, te)
lin_popt_vp, lin_pcov_vp = curve_fit(lin_fit, r, vp)

fig, axs = plt.subplots(1, 3, figsize=(8, 5))
axs = axs.flatten()

axs[0].scatter(r, ne)
axs[1].scatter(r, te)
axs[2].scatter(r, vp)

axs[0].plot(r, lin_fit(r, *lin_popt_ne))
axs[1].plot(r, lin_fit(r, *lin_popt_te))
axs[2].plot(r, lin_fit(r, *lin_popt_vp))

axs[0].plot(r, exp_fit(r - 2.32, *exp_popt_ne))
axs[1].plot(r, exp_fit(r - 2.32, *exp_popt_te))
axs[2].plot(r, exp_fit(r - 2.32, *exp_popt_vp))

axs[0].set_ylabel("ne (1e18 m-3)")
axs[1].set_ylabel("Te (eV)")
axs[2].set_ylabel("Vp (V)")

fig.tight_layout()
fig.show()
