import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import pretty_plots as pp


# Location of Excel file with "TOTAL DEPOSITION" table copy/pasted into it.
filename = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/test28_totdep.xlsx'

# Load the data into a Dataframe, drop the last row that has junk in it, and
# then rename the columns that is a string '0.1' to a float 0.0 (not sure why
# this happens).
print("Loading file...")
df = pd.read_excel(filename, index_col=0)
df.drop(df.columns[-1], axis=1, inplace=True)
df.rename({'0.1':0.0, '0.2':0.0}, axis=1, inplace=True)

# Get each side of probe, and a centerline (P=0).
itf_x = df.iloc[df.index > 0.0][0].index.values * 100
itf_y = df.iloc[df.index > 0.0][0].values * -1
otf_x = df.iloc[df.index < 0.0][0].index.values * -1 * 100
otf_y = df.iloc[df.index < 0.0][0].values * -1

def exp_fit(x, a, b):
    return a * np.exp(-x*b)

itf_popt, itf_pcov = curve_fit(exp_fit, itf_x, itf_y)
otf_popt, otf_pcov = curve_fit(exp_fit, otf_x, otf_y)
itf_x_fit = np.linspace(itf_x.min(), itf_x.max(), 100)
otf_x_fit = np.linspace(otf_x.min(), otf_x.max(), 100)
itf_y_fit = exp_fit(itf_x_fit, *itf_popt)
otf_y_fit = exp_fit(otf_x_fit, *otf_popt)

print('ITF Lambda = {:.2f}'.format(1/itf_popt[1]))
print('OTF Lambda = {:.2f}'.format(1/otf_popt[1]))

fig = None
fig = pp.pplot(itf_x_fit, itf_y_fit, '--', label='ITF Fit', color=6, fig=fig)
fig = pp.pplot(otf_x_fit, otf_y_fit, '--', label='OTF Fit', color=8, fig=fig)
fig = pp.pplot(itf_x, itf_y, '.', label='ITF', fig =fig)
fig = pp.pplot(otf_x, otf_y, '.', label='OTF', color=8, fig=fig)
