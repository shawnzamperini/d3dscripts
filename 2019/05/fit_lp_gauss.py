import pandas as pd
import numpy as np
import pretty_plots as pp
from scipy.special import erfc
from scipy.optimize import curve_fit


filename = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/outer_targ_data.xlsx'
df = pd.read_excel(filename, sheet_name='Sheet1')

fx = 5
psin = df['Psin'].values
ne   = df['ne (cm-3)'].values * 10**(-12)
te   = df['Te (eV)'].values * 10**(-1)
jsat = df['jsat (A/m2)'].values * 10**(-4)

# The ne fit's fine, so don't worry about it. For Te, lets just fit the peak
# to the convoluted Gaussian, but fit the parts outside to exponentials.
te_left = psin < 1.0; te_right = psin > 1.02
te_left_exp  = te[te_left]
te_right_exp = te[te_right]
peak = np.logical_and(psin > 1.0, psin < 1.02)
peak[np.where(peak == True)[0][0] - 3] = True

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

def gauss_conv_exp_fit(s, width, lambda_n, n0, n_bg, s0):
    return n0 / 2.0 * np.exp((width/(2*lambda_n*fx))**2 - (s-s0)/(lambda_n * fx)) * erfc(width/(2*lambda_n*fx) - (s-s0)/width)

popt_te_left,  pcov = curve_fit(exp_fit, psin[te_left],  te_left_exp,  maxfev=5000)
popt_te_right, pcov = curve_fit(exp_fit, psin[te_right], te_right_exp, maxfev=5000)

guess_ne = (0.05, 0.02, ne.max(), ne.min(), 1.0)
guess_te = (0.05, 0.02, te.max(), 0.0, 1.0)
guess_jsat = (0.05, 0.02, jsat.max(), 0.0, 1.0)
popt_ne, pcov_ne = curve_fit(gauss_conv_exp_fit, psin, ne, p0=guess_ne, maxfev=5000)
popt_te, pcov_te = curve_fit(gauss_conv_exp_fit, psin[peak], te[peak], p0=guess_te, maxfev=5000)
#popt_js, pcov_js = curve_fit(gauss_conv_exp_fit, psin, jsat, p0=guess_jsat, maxfev=5000)

ext_psin = np.linspace(0.96, 1.5, 100)
fig = pp.pplot(ext_psin, gauss_conv_exp_fit(ext_psin, *popt_ne), fmt='-')
fig = pp.pplot(psin, ne, fmt='.', xlabel='Psin', ylabel='ne', fig=fig)
fig = pp.pplot(psin, te, fmt='.')
fig = pp.pplot(psin[peak], gauss_conv_exp_fit(psin[peak], *popt_te), fmt='-', fig=fig)
fig = pp.pplot(psin[te_left], exp_fit(psin[te_left], *popt_te_left), '-', fig=fig)
fig = pp.pplot(psin[te_right], exp_fit(psin[te_right], *popt_te_right), '-', fig=fig)
#fig = pp.pplot(ext_psin, gauss_conv_exp_fit(ext_psin, *popt_js), fmt='-')
#fig = pp.pplot(psin, jsat, fmt='.', xlabel='Psin', ylabel='Jsat (A/m2)', fig=fig)

ext_x = np.linspace(psin.max(), 1.5, 20)

psin_out = np.append(psin[te_left], psin[peak])
psin_out = np.append(psin_out, psin[te_right])
psin_out = np.append(psin_out, ext_x)

left_out = exp_fit(psin[te_left], *popt_te_left)
peak_out = gauss_conv_exp_fit(psin[peak], *popt_te)
right_out = exp_fit(psin[te_right], *popt_te_right)
ext_out = exp_fit(ext_x, *popt_te_right)

te_out = np.append(left_out, peak_out)
te_out = np.append(te_out, right_out)
te_out = np.append(te_out, ext_out)

ne_out = gauss_conv_exp_fit(psin_out, *popt_ne)
