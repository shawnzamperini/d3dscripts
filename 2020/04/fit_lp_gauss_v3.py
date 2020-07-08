import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.optimize import curve_fit


# Fix the step in the LP profiles that occur for the floor probes right next to
# the shelf.
shelf_fix = False

# Multiply the right exponential value by a constant. Will automatically force
# the profiles to intersect at the right point.
te_right_mult = 0.5
jsat_right_mult = 3

# The Excel file with the actual LP data (outer target).
filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2020/04/lp_input_167247.xlsx'
df = pd.read_excel(filename, sheet_name='167247 LP')

# Values for fitting. Little multipliers sometimes help the fit, just don't
# forget to divide by them at the end.
te_mod = 10
jsat_mod = 10
psin = df['Psin'].values
te   = df['Te (eV)'].values / te_mod
jsat = df['jsat (A/cm2)'].values / jsat_mod

# Get rid of nans.
nans = np.isnan(te)
te   = te[~nans]
psin = psin[~nans]
jsat = jsat[~nans]

# Psin values determining the regions for the fits. Anything less than left or
# more than right are exponentially fit. In between them are fit to a
# convoluted Gaussian.
te_left_psin   = 1.0; te_right_psin   = 1.04
jsat_left_psin = 1.004; jsat_right_psin = 1.04
te_left   = psin < te_left_psin;   te_right   = psin > te_right_psin
jsat_left = psin < jsat_left_psin; jsat_right = psin > jsat_right_psin
te_peak   = np.logical_and(psin >= te_left_psin,   psin <= te_right_psin)
jsat_peak = np.logical_and(psin >= jsat_left_psin, psin <= jsat_right_psin)

# The fitting functions.
def exp_fit(x, a, b):
    return a * np.exp(-b * x)

def gauss_conv_exp_fit(s, width, lambda_n, n0, n_bg, s0):
    fx = 5
    return n0 / 2.0 * np.exp((width/(2*lambda_n*fx))**2 - (s-s0)/(lambda_n *
      fx)) * erfc(width/(2*lambda_n*fx) - (s-s0)/width)

# Exponential fits.
popt_te_left, pcov    = curve_fit(exp_fit, psin[te_left],    te[te_left],      p0=(1, 10), maxfev=50000)
popt_te_right, pcov   = curve_fit(exp_fit, psin[te_right],   te[te_right],     p0=(1, 10), maxfev=50000)
popt_jsat_left, pcov  = curve_fit(exp_fit, psin[jsat_left],  jsat[jsat_left],  p0=(1, 10), maxfev=50000)
popt_jsat_right, pcov = curve_fit(exp_fit, psin[jsat_right], jsat[jsat_right], p0=(1, 10), maxfev=50000)

# Convoluted Gaussian fit to Te.
guess_te = (0.05, 0.02, te.max(), 0.0, 1.0)
popt_te, pcov_te = curve_fit(gauss_conv_exp_fit, psin[te_peak], te[te_peak],
                     p0=guess_te, maxfev=50000)

# Convoluted Gaussian fit to jsat.
guess_jsat = (0.05, 0.02, jsat.max(), 0.0, 1.0)
popt_jsat, pcov_jsat = curve_fit(gauss_conv_exp_fit, psin[jsat_peak],
                         jsat[jsat_peak], p0=guess_jsat, maxfev=5000)

# Will combine the fits as 100 data points.
psin_fit = np.linspace(psin.min(), psin.max(), 150)
jsat_fit = np.zeros(len(psin_fit))
te_fit   = np.zeros(len(psin_fit))

# Fill in the left and right exponential fit data.
te_fit[psin_fit < te_left_psin]      = exp_fit(psin_fit[psin_fit < te_left_psin],    *popt_te_left)
te_fit[psin_fit > te_right_psin]     = exp_fit(psin_fit[psin_fit > te_right_psin],   *popt_te_right)
jsat_fit[psin_fit < jsat_left_psin]  = exp_fit(psin_fit[psin_fit < jsat_left_psin],  *popt_jsat_left)
jsat_fit[psin_fit > jsat_right_psin] = exp_fit(psin_fit[psin_fit > jsat_right_psin], *popt_jsat_right)

# Fill in the Gaussian data.
peak_te_fit = np.logical_and(psin_fit>=te_left_psin, psin_fit<=te_right_psin)
te_fit[peak_te_fit] = gauss_conv_exp_fit(psin_fit[peak_te_fit], *popt_te)
peak_jsat_fit = np.logical_and(psin_fit>=jsat_left_psin, psin_fit<=jsat_right_psin)
jsat_fit[peak_jsat_fit] = gauss_conv_exp_fit(psin_fit[peak_jsat_fit], *popt_jsat)

# Apply shelf fix procedure to the Te data.
if shelf_fix:
    te_fit_org = np.copy(te_fit)
    # Extend the right exponential to the left until it's less than the Gaussian
    # fit (i.e. it intersects it).
    for i in range(0, len(psin_fit)):
        new_te = exp_fit(psin_fit[-(i+1)], *popt_te_right)
        if new_te < gauss_conv_exp_fit(psin_fit[-(i+1)], *popt_te):
            break
        te_fit[-(i+1)] = new_te

# Multiply the right half of the exponential data by a constant, making sure the
# data intersects at the right point.
if te_right_mult is not None:
    te_fit_old = np.copy(te_fit)
    for i in range(0, len(psin_fit)):
        new_te = exp_fit(psin_fit[-(i+1)], *popt_te_right) * te_right_mult
        if new_te < gauss_conv_exp_fit(psin_fit[-(i+1)], *popt_te):
            break
        te_fit[-(i+1)] = new_te

# Same for Jsat.
if jsat_right_mult is not None:
    jsat_fit_old = np.copy(jsat_fit)
    for i in range(0, len(psin_fit)):
        new_jsat = exp_fit(psin_fit[-(i+1)], *popt_jsat_right) * jsat_right_mult
        if new_jsat < gauss_conv_exp_fit(psin_fit[-(i+1)], *popt_jsat):
            break
        jsat_fit[-(i+1)] = new_jsat

# End od the day multiply by the number we divided by.
te_fit *= te_mod
jsat_fit *= jsat_mod

# Text for plots.
jsat_str = 'Left:  {:4.2e} * exp({:.2e}*x)\n'.format(popt_jsat_left[0],  -popt_jsat_left[1]) + \
           'Right: {:4.2e} * exp({:.2e}*x)'.format(popt_jsat_right[0], -popt_jsat_right[1])
te_str   = 'Left:  {:4.2e} * exp({:.2e}*x)\n'.format(popt_te_left[0],  -popt_te_left[1]) + \
           'Right: {:4.2e} * exp({:.2e}*x)'.format(popt_te_right[0], -popt_te_right[1])

# Plots.
fig, axs = plt.subplots(1, 2, figsize=(15, 8))
ax1 = axs.flatten()[0]; ax2 = axs.flatten()[1]
ax1.plot(psin, te * te_mod, 'k.')
ax1.plot(psin_fit, te_fit, 'r-')
if shelf_fix:
    ax1.plot(psin_fit, te_fit_org*te_mod, 'r--')
ax1.set_xlabel('Psin')
ax1.set_ylabel('Te (eV)')
ax1.annotate(te_str, xy=(0.5, 0.75), xycoords='axes fraction')
ax2.plot(psin, jsat * jsat_mod, 'k.')
ax2.plot(psin_fit, jsat_fit, 'b-')
ax2.set_xlabel('Psin')
ax2.set_ylabel('Jsat (A/cm2)')
ax2.annotate(jsat_str, xy=(0.5, 0.75), xycoords='axes fraction')
fig.tight_layout()
fig.show()
