import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


filename = '/home/shawn/d3dscripts/Data/LModeProbes.xlsx'
df       = pd.read_excel(filename, sheet_name='A2', usecols=range(0,14)).dropna()

xd     = df['rminrsep_D']
xd_err = df['rminrsep_err_D']
yd     = df['w_areal_D']
yd_err = df['w_areal_err_D']
xu     = df['rminrsep_U']
xu_err = df['rminrsep_err_U']
yu     = df['w_areal_U']
yu_err = df['w_areal_err_U']

def exp_fit(x, a, b):
    return a * np.exp(-x*b)

poptd, pcovd = curve_fit(exp_fit, xd[:-3], yd[:-3])
x_fitd = np.linspace(7, 15, 100)
y_fitd = exp_fit(x_fitd, *poptd)
poptu, pcovu = curve_fit(exp_fit, xu[:-3], yu[:-3])
x_fitu = np.linspace(7, 15, 100)
y_fitu = exp_fit(x_fitu, *poptu)

font = {'fontsize' : 24, 'weight' : 'bold'}
plt.style.use('seaborn')
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.errorbar(xd, yd, yd_err, xd_err, 'r.', ms=20, capsize=5, capthick=1, label='A2-ITF')
ax1.errorbar(xu[:-1], yu[:-1], yu_err[:-1], xu_err[:-1], 'b.', ms=20, capsize=5, capthick=1, label='A2-OTF')
ax1.plot(x_fitd, y_fitd, 'r--')
ax1.plot(x_fitu, y_fitu, 'b--')
ax1.set_xlabel(r'$\mathrm{\bf{R-R_{sep}\ (cm)}}$', font)
ax1.set_ylabel(r'$\mathrm{\bf{W\ Density\ (10^{18}\ cm^{-2})}}$', font)
ax1.legend(frameon=True, fontsize=24)
ax1.set_ylim([None,0.6])
ax1.tick_params(labelsize=22)
fig.tight_layout()
fig.show()
