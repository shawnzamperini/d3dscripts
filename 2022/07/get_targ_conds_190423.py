# Script to output per-ring specification for the outer target conditions
# for the 190423 grid.
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import oedge_plots
from scipy.optimize import curve_fit


# Load the LP data.
xlpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/lp_190423.xlsx"
xl = pd.read_excel(xlpath, sheet_name="lp")[["fit_psin", "fit_te", "fit_jsat"]].dropna()
lp_psin = xl["fit_psin"].values
lp_te = xl["fit_te"].values
lp_jsat = xl["fit_jsat"].values

# Create interpolation from LP data, anything outside of the bounds just use
# the last value.
f_te = interp1d(lp_psin, lp_te, bounds_error=False, fill_value=(lp_te[0], lp_te[-1]))
f_jsat1 = interp1d(lp_psin, lp_jsat, bounds_error=False, fill_value=(lp_jsat[0], lp_jsat[-1]))
f_jsat2 = interp1d(lp_psin, lp_jsat, bounds_error=True)

# If outside of the interpolation, then just use an exponential fit to
# extrapolate out for jsat only (Te seems to be fine staying the same).
def exp_fit(x, a, b):
    return a * np.exp(-b * x)
exp_fit_psin = np.linspace(1.01, 1.02, 50)
popt, pcov = curve_fit(exp_fit, exp_fit_psin-1, f_jsat2(exp_fit_psin)/1e5)
int_psin = np.linspace(1.02, 1.10, 100) - 1
f_jsat_exp = interp1d(int_psin + 1, exp_fit(int_psin, *popt)*1e5)


# Load the grid data, in particular the rings and psin of each ring.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d_extgridtest.nc"
op = oedge_plots.OedgePlots(ncpath)
nrs = int(op.nc["NRS"][:])
rings = np.arange(1, nrs+1, dtype=int)
ring_psins = []; ring_tes = []; ring_jsats = []
for i in range(0, len(rings)):

    # Make sure to zero index.
    psin = float(op.nc["PSIFL"][rings[i]-1][0])
    ring_psins.append(psin)

    # Get corresponding Te, jsat data for the rings.
    ring_tes.append(f_te(psin))

    # If in the PFZ just extrapolate, don't really care here since we'll manually change it.
    if psin < 1:
        ring_jsats.append(f_jsat1(psin))
    else:
        try:
            ring_jsats.append(f_jsat2(psin))
        except:
            ring_jsats.append(f_jsat_exp(psin))

# Print out in a format for the input file.
for i in range(0, len(rings)):
    print("{:4d} {:5.2f} {:5.2f} {:5.2e}".format(rings[i], ring_tes[i], ring_tes[i], ring_jsats[i]))
