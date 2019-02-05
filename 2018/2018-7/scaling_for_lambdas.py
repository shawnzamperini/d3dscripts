from scipy.optimize import curve_fit
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt


filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018-7/scaling.xlsx'
df = pd.read_excel(filename, sheet_name='Lambda Trends')

# Only use H-mode shots.
df.drop([0,1,2,3,4], inplace=True)

lamb_rat = df['Lambda Ratio']
fg = df['Greenwald Fraction']
q95 = df['q95']

def fit_func(X1X2, a, b, c):
    x1, x2 = X1X2
    return a * x1**b * x2**c

popt, pcov = curve_fit(fit_func, (fg, q95), lamb_rat)
fit_vals = fit_func((fg, q95), *popt)

textstr = r'$\mathrm{R_{scaling} = a\ *\ f_G^{b}\ *\ Q95^{c}}$' + \
          '\na = {:.2f}'.format(popt[0]) + \
          '\nb = {:.2f}'.format(popt[1]) + \
          '\nc = {:.2f}'.format(popt[2])

# Some properties for the plots.
props = dict(alpha=0.25, facecolor='k')
font = {'fontsize':24, 'weight':'bold'}
plt.style.use('seaborn')
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(fit_vals, lamb_rat, 'k.', ms=10)
ax1.set_xlabel(r'$\mathrm{R_{scaling}}$', font)
ax1.set_ylabel(r'$\mathrm{R_{measured}}$', font)
ax1.plot(range(0,5), 'k--')
ax1.set_xlim([0.25,1])
ax1.set_ylim([0.25,1])
ax1.set_title('ITF/OTF Lambda Ratio', font)
plt.tick_params(axis='both', which='major', labelsize=18)
ax1.text(0.6, 0.1, textstr, transform=ax1.transAxes, fontsize=18, bbox=props)
fig.tight_layout()
fig.show()
