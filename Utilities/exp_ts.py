import pretty_plots as pp
import numpy as np
from ThomsonClass import ThomsonClass
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from gadata import gadata
import MDSplus as mds


# Enable interactive mode so user can input with graphs still up.
plt.ion()

def ts_fitting(shot, tmin, tmax, tmany, tree):

    # Load the TS data.
    ts = ThomsonClass(shot, 'core')
    #ts = ThomsonClass(shot, 'divertor')
    ts.load_ts()
    ts.map_to_efit(np.linspace(tmin, tmax, tmany), tree=tree, trunc_div=False)

    # Pull out the arrays.
    r  = ts.avg_omp['RminRsep_omp'] * 100
    te = ts.avg_omp['Te_omp']
    ne = ts.avg_omp['ne_omp'] * 10**(-18)
    r_err  = ts.avg_omp['RminRsep_omp_err'] * 100  # m to cm
    te_err = ts.avg_omp['Te_omp_err']
    ne_err = ts.avg_omp['ne_omp_err'] * 10**(-18) # m-3 to 10^-18 m-3

    # Plot with chord numbers next to each so you can point out the bad ones.
    fig = plt.figure(figsize=(10,5))
    ax_te = fig.add_subplot(121)
    ax_ne = fig.add_subplot(122)
    ax_te.errorbar(r, te, te_err, r_err, fmt='k.', ms=8)
    ax_ne.errorbar(r, ne, ne_err, r_err, fmt='k.', ms=8)
    ax_te.set_xlabel('R-Rsep OMP(cm)')
    ax_ne.set_xlabel('R-Rsep OMP(cm)')
    ax_te.set_ylabel('Te (eV)')
    ax_ne.set_ylabel(r'$\mathrm{ne\ (10^{18}\ m^{-3})}$')
    for i, chord in enumerate(np.arange(0, len(r))):
        ax_te.annotate(chord, (r[i], te[i]))
        ax_ne.annotate(chord, (r[i], ne[i]))
    fig.tight_layout()
    fig.show()

    while True:

        # Choose chords for fitting. Just start with the first three inside the sep
        # and the rest outwards.
        num_in_include = 3
        first_out = np.where(r > 0)[0][-1]
        include = num_in_include + first_out

        # Zoom into the fitting region on the plots.
        ax_te.set_xlim(1.1 * r[:include].min(), 1.1 * r[:include].max())
        ax_te.set_ylim(1.1 * te[:include].min(), 1.1 * te[:include].max())
        ax_ne.set_xlim(1.1 * r[:include].min(), 1.1 * r[:include].max())
        ax_ne.set_ylim(1.1 * ne[:include].min(), 1.1 * ne[:include].max())

        # Ask if any extra chords should be excluded.
        print('\nChords for fitting: ', end=''); print(*np.arange(0, len(r))[:include])
        exclude = input('Chords to exclude (separated by commas, press enter if none): ').split(',')

        if exclude != ['']:
            r_tofit  = np.delete(r[:include],  exclude)
            te_tofit = np.delete(te[:include], exclude)
            ne_tofit = np.delete(ne[:include], exclude)
        else:
            r_tofit  = np.array(r[:include])
            te_tofit = np.array(te[:include])
            ne_tofit = np.array(ne[:include])

        def exp_fit(x, a, b):
            return a * np.exp(-x*b)

        # Perform the fitting.
        popt_te, pcov_te = curve_fit(exp_fit, r_tofit, te_tofit, maxfev=10000)
        popt_ne, pcov_ne = curve_fit(exp_fit, r_tofit, ne_tofit, maxfev=10000)

        r_fit = np.linspace(r_tofit.min(), r_tofit.max(), 50)
        te_fit = exp_fit(r_fit, *popt_te)
        ne_fit = exp_fit(r_fit, *popt_ne)

        # Plot the fit.
        ax_te.plot(r_fit, te_fit, 'k--', lw=5)
        ax_ne.plot(r_fit, ne_fit, 'k--', lw=5)
        fig.show()

        print('Te: {:.2f} * exp(-r / {:.2f})'.format(popt_te[0], 1/popt_te[1]))
        print('ne: {:.2f} * exp(-r / {:.2f})'.format(popt_ne[0], 1/popt_ne[1]))

        print('\nTe falloff (cm): {:.2f}'.format(1 / popt_te[1]))
        print('ne falloff (cm): {:.2f}'.format(1 / popt_ne[1]))

        ans = input('Try fitting again? (y/n): ')
        if ans == 'n':
            break

    return r, te, ne

def flat_top(shot):

    # Load densv2 tag data.
    conn = mds.Connection('localhost')
    ga_obj = gadata('densv2', shot, connection=conn)
    time = ga_obj.xdata
    dens = ga_obj.zdata

    # Plot and ask for the flat top range.
    fig = pp.pplot(time, dens, fmt='-',  xlabel='Time (ms)', ylabel=r'$\mathrm{\bar{n_e}\ (m^{-3})}$')
    minmax = input('Enter time min/max for analysis range (separated by commas): ').split(',')

    # Return requested time range for use in TS function (min, max).
    return int(minmax[0]), int(minmax[1])

def run(shot, tmany=10, tree='EFIT04'):

    min, max = flat_top(shot)
    r, te, ne = ts_fitting(shot, min, max, tmany, tree=tree)

    return {'r':r, 'te':te, 'ne':ne}
