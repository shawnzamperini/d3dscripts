from scipy.optimize import curve_fit
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt


filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018-5/max_w_vs_flux.xlsx'
df = pd.read_excel(filename, sheet_name='Total Content', skiprows=[0])

df = df.iloc[:13]

asymm     = np.array(df['Ratio'].values, dtype=np.float64)[4:15]
asymm_err = np.array(df['Error'].values, dtype=np.float64)[4:15]
lambda_ne = np.array(df['Density Fall Off (cm)'].values, dtype=np.float64)[4:15]
ng        = np.array(df['Greenwald Fraction'].values, dtype=np.float64)[4:15]
ip        = np.array(df['IP (MA)'].values, dtype=np.float64)[4:15]
ne_sep    = np.array(df['ne sep (m-3)'].values, dtype=np.float64)[4:15]
balloon   = np.array(df['Ballooning Parameter'].values, dtype=np.float64)[4:15]

# Rescale all between 0 and 1.
#asymm = np.interp(asymm, (asymm.min(), asymm.max()), (0,1))
lambda_ne = np.interp(lambda_ne, (lambda_ne.min(), lambda_ne.max()), (0,1))
ng = np.interp(ng, (ng.min(), ng.max()), (0,1))

"""
asymm     = np.array(df['Ratio'].values, dtype=np.float64)
asymm_err = np.array(df['Error'].values, dtype=np.float64)
lambda_ne = np.array(df['Density Fall Off (m)'].values, dtype=np.float64)
ng        = np.array(df['Greenwald Fraction'].values, dtype=np.float64)
ip        = np.array(df['IP (MA)'].values, dtype=np.float64)
ne_sep    = np.array(df['ne sep (1e18 m-3)'].values, dtype=np.float64)
"""
def fit_func(X1X2, a, b, c):
    x1, x2 = X1X2
    return c * x1**a * x2**b

popt, pcov = curve_fit(fit_func, (lambda_ne, ng), asymm)

fit_vals = fit_func((lambda_ne, ng), *popt)

textstr = r'$\mathrm{R_{fit} = A\ *\ \lambda_{ne}^{B}\ *\ f_g^{C}}$' + \
          '\nA = {:.2f}'.format(popt[2]) + \
          '\nB = {:.2f}'.format(popt[0]) + \
          '\nC = {:.2f}'.format(popt[1])
#textstr = r'$\mathrm{fit = c\ *\ \lambda_{ne}^{a}\ *\ n_g^{b}\na={:.2f}\nb={:.2f}\nc={:.2f}}$'.format(popt[1], popt[2], popt[0])
props = dict(alpha=0.25, facecolor='k')

font = {'fontsize' : 24,
        'weight'   : 'bold'}
plt.style.use('seaborn')

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(x=fit_vals, y=asymm, yerr=asymm_err, ms=10, fmt='k.', capthick=1, capsize=3)
ax1.set_xlabel(r'$\mathrm{R_{fit}}$', font)
ax1.set_ylabel(r'$\mathrm{R_{measured}}$', font)
ax1.plot(range(0,5), 'k--')
ax1.set_xlim([0,2])
ax1.set_ylim([0,2])
ax1.set_title('ITF/OTF Total Content Ratio', font)
plt.tick_params(axis='both', which='major', labelsize=18)
ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=18, bbox=props)
fig.tight_layout()
fig.show()


"""
# Try the scaling with lambda_ne and Ip now (Ip is simpler).
popt, pcov = curve_fit(fit_func, (lambda_ne, ip), asymm)
fit_vals = fit_func((lambda_ne, ip), *popt)
textstr = r'$\mathrm{fit = c\ *\ \lambda_{ne}^{a}\ *\ Ip^{b}}$' + \
          '\na = {:.2f}'.format(popt[0]) + \
          '\nb = {:.2f}'.format(popt[1]) + \
          '\nc = {:.2f}'.format(popt[2])

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(x=fit_vals, y=asymm, yerr=asymm_err, ms=10, fmt='k.', capthick=1, capsize=3)
ax1.set_xlabel('Fit')
ax1.set_ylabel('Experimental')
ax1.plot(range(0,5), 'k--')
ax1.set_xlim([0,4])
ax1.set_ylim([0,4])
ax1.set_title('ITF/OTF Total Content Ratio')
ax1.text(0.7, 0.1, textstr, transform=ax1.transAxes, fontsize=12, bbox=props)
fig.tight_layout()
fig.show()

# Try the scaling with lambda_ne and balloon now.
popt, pcov = curve_fit(fit_func, (lambda_ne, balloon), asymm)
fit_vals = fit_func((lambda_ne, balloon), *popt)
textstr = r'$\mathrm{fit = c\ *\ \lambda_{ne}^{a}\ *\ \alpha^{b}}$' + \
          '\na = {:.2f}'.format(popt[0]) + \
          '\nb = {:.2f}'.format(popt[1]) + \
          '\nc = {:.2f}'.format(popt[2])

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(x=fit_vals, y=asymm, yerr=asymm_err, ms=10, fmt='k.', capthick=1, capsize=3)
ax1.set_xlabel('Fit')
ax1.set_ylabel('Experimental')
ax1.plot(range(0,5), 'k--')
ax1.set_xlim([0,4])
ax1.set_ylim([0,4])
ax1.set_title('ITF/OTF Total Content Ratio')
ax1.text(0.7, 0.1, textstr, transform=ax1.transAxes, fontsize=12, bbox=props)
fig.tight_layout()
fig.show()

# Do all 4 just for fun.
def fit_func(X1X2X3X4, a, b, c, d, e):
    x1, x2, x3, x4 = X1X2X3X4
    return e * x1**a * x2**b * x3**c * x4**d

popt, pcov = curve_fit(fit_func, (lambda_ne, ng, ip, ne_sep), asymm)
fit_vals = fit_func((lambda_ne, ng, ip, ne_sep), *popt)
textstr = r'$\mathrm{fit = e\ *\ \lambda_{ne}^{a}\ *\ fg^{b}}\ *\ Ip^c\ *\ ne_{sep}^d$' + \
          '\na = {:.2f}'.format(popt[0]) + \
          '\nb = {:.2f}'.format(popt[1]) + \
          '\nc = {:.2f}'.format(popt[2]) + \
          '\nd = {:.2f}'.format(popt[3]) + \
          '\ne = {:.2f}'.format(popt[4])

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(x=fit_vals, y=asymm, yerr=asymm_err, ms=10, fmt='k.', capthick=1, capsize=3)
ax1.set_xlabel('Fit')
ax1.set_ylabel('Experimental')
ax1.plot(range(0,5), 'k--')
ax1.set_xlim([0,4])
ax1.set_ylim([0,4])
ax1.set_title('ITF/OTF Total Content Ratio')
ax1.text(0.7, 0.1, textstr, transform=ax1.transAxes, fontsize=12, bbox=props)
fig.tight_layout()
fig.show()
"""
