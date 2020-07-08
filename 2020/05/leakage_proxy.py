import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Some constants.
lp_shift  = -0.003  # Shift LP data in psin.
fit_shift = -0.003  # Shift fit data in psin.

# Load data from the specific sheet.
file_path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Slides, Sheets and Documents/2020/04/lp_input_167247.xlsx'
df = pd.read_excel(file_path, sheet_name='167247 Python')

# Helper function to get rid of nans.
def get_array(col_name):
    vals = df[col_name].values
    vals = vals[~np.isnan(vals)]
    return vals

# Individual arrays for everything.
lp_psin  = get_array('LP psin') + lp_shift
lp_te    = get_array('LP Te')
fit_psin = get_array('Fit psin') + fit_shift
fit_te   = get_array('Fit Te')
ts_psin  = get_array('TS psin')
ts_te    = get_array('TS Te')
ts1_psin = get_array('TS psin.1')  # 167277
ts1_te   = get_array('TS Te.1')    # 167277

# Interpolation functions so everything can use the same psins.
f_lp  = interp1d(lp_psin, lp_te, 'cubic')
f_fit = interp1d(fit_psin, fit_te, 'cubic')
f_ts  = interp1d(ts_psin, ts_te, 'cubic')
f_ts1 = interp1d(ts1_psin, ts1_te, 'cubic')

# A common set of psin values to use.
low  = max(lp_psin.min(), fit_psin.min(), ts_psin.min())
high = min(lp_psin.max(), fit_psin.max(), ts_psin.max())
#psin = np.linspace(1.0, 1.1, 50)
psin = np.linspace(low, high, 500)

# Use interpolation to get values all at the same psins.
lp  = f_lp(psin)
fit = f_fit(psin)
ts  = f_ts(psin)
ts1 = f_ts1(psin)

# Estimate the leakage proxy.
leak_lp  = (ts - lp)
leak_fit = (ts - fit)
leak_lp1 = (ts1 - lp)
leak_fit1 = (ts1 - fit)

# Plot.
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, figsize=(15, 5))
ax1.plot(psin, leak_lp, label='167247')
ax2.plot(psin, leak_fit, label='167247')
ax1.plot(psin, leak_lp1, label='167277')
ax2.plot(psin, leak_fit1, label='167277')
ax3.plot(psin, leak_fit1/leak_fit)
ax1.set_xlabel('Psin')
ax2.set_xlabel('Psin')
ax1.set_ylabel('Leakage Proxy (Tu - Td)')
ax3.set_ylabel('Proxy 167277/167247')
#ax1.set_ylim([0, 3.5])
#ax2.set_ylim([0, 3.5])
ax3.set_ylim([0, 6])
ax3.set_xlim([None, 1.08])
ax1.legend()
ax2.legend()
fig.tight_layout()
fig.show()
