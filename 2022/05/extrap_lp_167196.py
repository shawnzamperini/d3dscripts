# From the Lanmguir probe data, extrapolate it out past where it ends with an exponential fit.
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



datapath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/data_for_167196.xlsx"
lp = pd.read_excel(datapath, sheet_name="lp_data_v2").replace("", np.nan)

# This values have already been filtered a bit, so we fit to them.
lp_psin = lp["psin.1"].values
lp_ne = lp["ne"].values  # not a typo

# The exponential region is only the trailing points on the outside.
mask = np.logical_and(lp_psin > 1.04, lp_psin < 1.09)
psin_fit = lp_psin[mask]
ne_fit = lp_ne[mask]

def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c  # We manually set the offset, arbitrarily mostly just so it doesn't go negative.

# Offset the data so it starts at zero, make ne smaller, good for fit.
psin_off = psin_fit - psin_fit.min()
ne_off = ne_fit / 1e18
popt_ne, pcov_ne = curve_fit(exp_fit, psin_off, ne_off)

# See what the fit predicts past the target LP data.
mask = lp_psin > 1.04
psin_extrap = lp_psin[mask] - psin_fit.min()
ne_extrap = exp_fit(psin_extrap, *popt_ne) * 1e18
psin_extrap_plot = psin_extrap + psin_fit.min()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), sharex=True)
ax1.scatter(lp_psin, lp_ne)
ax1.scatter(psin_fit, ne_fit)
ax1.scatter(psin_extrap_plot, ne_extrap)
ax1.set_xlim([0.8, 1.36])
fig.tight_layout()
fig.show()

# Print out the data so it can be copy/pasted into the Excel sheet.
print("Psin values")
for i in psin_extrap_plot:
    print(i)
print("\nne values")
for i in ne_extrap:
    print(i)
