# This script provides a plot of input Te data to 3DLIM to demonstrate how we
# come up with it.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/ts167247_final_v2.xlsx"
ts_df = pd.read_excel(ts_path)

ts_core = ts_df[ts_df["System"] == "core"]
psin = ts_core["Psin"]
te = ts_core["Te (eV)"]
ne = ts_core["ne (m-3)"] * 1e-18

mask = psin >= 1.0
psin = psin[mask]
te = te[mask]
ne = ne[mask]

#mask = np.logical_and(te > 1, te < 80)
#psin = psin[mask]
#te = te[mask]

# Convert to R-Rsep OMP in cm.
rmrsomp = 42.239 * psin - 43.265
rmrsomp  = rmrsomp - rmrsomp.min()

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

popt, pcov = curve_fit(exp_fit, rmrsomp, ne)

x = np.arange(0, 18)
fit_ne = exp_fit(x, *popt)
fit_ne[fit_ne < 2] = 2

fig, ax = plt.subplots()
ax.scatter(rmrsomp, ne)
ax.plot(x, fit_ne)
ax.set_yscale("log")
fig.tight_layout()
fig.show()
