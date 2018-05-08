"""
Plotting the flux at the tip of the probes with max W content.
"""

import get_ts            as ts
import numpy             as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if True:
    ts_dict = ts.run_script(167320, 'core', tmin=2000, tmax=4000)

x_ts = ts_dict['psins']['avg_omps'] * 100
y_te = ts_dict['psins']['avg_Tes']
y_ne = ts_dict['psins']['avg_nes'] * 10e-18

x_ts = x_ts[np.where(x_ts > 0)[0]]
y_te = y_te[np.where(x_ts > 0)[0]]
y_ne = y_ne[np.where(x_ts > 0)[0]]

m_deut  = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
y_cs    = np.sqrt(2*y_te/m_deut)
y_flux  = 0.5 * y_ne * 10e18 * y_cs
y_flux *= 10e-25

def exp_fit(x, a, b):
    return a * np.exp(-x*b)

popt_te,   pcov_te     = curve_fit(exp_fit, x_ts, y_te)
popt_ne,   pcov_ne     = curve_fit(exp_fit, x_ts, y_ne)
popt_flux, pcov_flux   = curve_fit(exp_fit, x_ts, y_flux)
x_fit            = np.linspace(0, 14, 100)
y_fit_te         = exp_fit(x_fit, *popt_te)
y_fit_ne         = exp_fit(x_fit, *popt_ne)
y_fit_flux       = exp_fit(x_fit, *popt_flux)

def plot_flux():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.semilogy(x_ts, y_flux*10e25, 'k.')
    ax1.semilogy(x_fit, y_fit_flux*10e25, 'k--')
    ax1.set_xlabel('R-Rsep omp (cm)')
    ax1.set_ylabel('Plasma Flux m-2 s-1')
    plt.show()

plot_flux()

def flux_at_tip(r_tip):
    return exp_fit(r_tip, *popt_flux) * 10e25
