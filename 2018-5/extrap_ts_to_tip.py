"""
Plotting the flux at the tip of the probes with max W content.
"""

import get_ts            as ts
import numpy             as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys


"""
lim value for each shot (some of them). Goal is to fit data outside the steep
drop region, but not so far out that the data become unreliable.

shot   : lowlim, uplim
167219 : 1, None
167277 : None, -6
167279 : None, -9
167320 : None, None
167321 : None, None
167322 : None, -5
167353 : None, None
167405 : None, -5
167408 : None, None
167463 : None, None
167481 : None, -5
167530 : 5, None
167531 : 3, None
167534 : None, None
167536 : None, -4


"""
lowlim = None     # Furthest points out to ignore. Must be positive.
uplim  = None  # Closest points to ignore. Must be negative.

if True:
    ts_dict = ts.run_script(167481, 'core', tmin=2000, tmax=4000)

x_ts = ts_dict['psins']['avg_omps'] * 100
y_te = ts_dict['psins']['avg_Tes']
y_ne = ts_dict['psins']['avg_nes'] * 10e-18

x_ts_err = ts_dict['psins']['avg_omps_err'] * 100
y_te_err = ts_dict['psins']['avg_Tes_err']
y_ne_err = ts_dict['psins']['avg_nes_err'] * 10e-18

x_ts = x_ts[np.where(x_ts > 0)[0]][lowlim:uplim]
y_te = y_te[np.where(x_ts > 0)[0]]
y_ne = y_ne[np.where(x_ts > 0)[0]]

x_ts_err = x_ts_err[np.where(x_ts > 0)[0]]
y_te_err = y_te_err[np.where(x_ts > 0)[0]]
y_ne_err = y_ne_err[np.where(x_ts > 0)[0]]

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
    ax1 = fig.add_subplot(211)
    ax1.plot(x_ts, y_flux*10e25, 'k.')
    ax1.plot(x_fit, y_fit_flux*10e25, 'k--')
    #ax1.set_xlabel('R-Rsep omp (cm)')
    ax1.set_ylabel(r'Plasma Flux $\mathrm{m^{-2} s^{-1}}$')
    ax2 = fig.add_subplot(212)
    ax2.plot(x_ts, y_ne, 'k.')
    ax2.plot(x_fit, y_fit_ne, 'k--')
    ax2.set_xlabel('R-Rsep omp (cm)')
    ax2.set_ylabel(r'Density (1e18 $\mathrm{m^{-3}}$)')
    fig.tight_layout()
    plt.show()

#plot_flux()

print("Density fall off length: {:.3f}".format(1/popt_ne[1]))
errs = np.sqrt(np.diag(pcov_ne))
print("Density fall off length error: {:.4f}".format((errs[1]/popt_ne[1]) * (1/popt_ne[1])))
ne_sep = exp_fit(0, *popt_ne)
print("Density at separatrix: {:.4f}".format(ne_sep))

print("Temperature fall off length: {:.3f}".format(1/popt_te[1]))
errs = np.sqrt(np.diag(pcov_te))
print("Temperature fall off length error: {:.4f}".format((errs[1]/popt_te[1]) * (1/popt_te[1])))
te_sep = exp_fit(0, *popt_te)
print("Temperature at separatrix: {:.4f}".format(te_sep))

def flux_at_tip(r_tip):
    return exp_fit(r_tip, *popt_flux) * 10e25

if False:
    font = {'fontsize' : 24,
            'weight'   : 'bold'}
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(x_ts, y_ne, 'k.', ms=20)
    ax1.plot(x_fit, y_fit_ne, 'k--', lw=3)
    ax1.set_xlabel('R-Rsep omp (cm)', font)
    ax1.set_ylabel(r'$\mathrm{\mathbf{Density\ (1e18\ m^{-3}}}$)', font)
    ax1.tick_params(labelsize=14)
    fig.tight_layout()
    plt.show()

if True:
    font = {'fontsize' : 24,
            'weight'   : 'bold'}
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.errorbar(x_ts, y_ne, y_ne_err, 0.5, 'k.', ms=20, capsize=5, capthick=1)
    ax1.plot(x_fit, y_fit_ne, 'k--', lw=3)
    ax1.set_xlabel('R-Rsep omp (cm)', font)
    ax1.set_ylabel(r'$\mathrm{\mathbf{Density\ (1e18\ m^{-3}}}$)', font)
    ax1.tick_params(labelsize=22)
    fig.tight_layout()
    plt.show()

if True:
    font = {'fontsize' : 24,
            'weight'   : 'bold'}
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.errorbar(x_ts, y_te, y_te_err, 0.5, 'k.', ms=20, capsize=5, capthick=1)
    ax1.plot(x_fit, y_fit_te, 'k--', lw=3)
    ax1.set_xlabel('R-Rsep omp (cm)', font)
    ax1.set_ylabel(r'$\mathrm{\mathbf{T_e (eV)}}$)', font)
    ax1.tick_params(labelsize=22)
    fig.tight_layout()
    plt.show()
