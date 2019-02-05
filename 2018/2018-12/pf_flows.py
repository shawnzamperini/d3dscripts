import get_lp as lp
import numpy as np
import pretty_plots as pp
from scipy.optimize import curve_fit
from gadata import gadata
import MDSplus
from ThomsonClass import ThomsonClass


def use_lp(shot, tmin, tmax, p_exclude=None, V_adjust=False):
    """
    Calculate the PF stuff using LP data (i.e. near the source). Typically
    less data points available and the integrity suffers. Especially when
    the SP is on the floor.

    p_exclude: Exclude these first (last) probes in plots using a positive
               (negative) value for an index.
    """

    # Get the LP data.
    lps = lp.get_dict_of_lps(shot)

    # Get the LP data into plottable arrays.
    Rs  = np.array([])
    Tes = np.array([])
    nes = np.array([])
    Vs  = np.array([])
    for lp_num in lps.keys():
        #print('lp: {}'.format(lp_num))

        # Don't want the upper divertor or anything.
        if lps[lp_num]['pnum'] < 39:
            tmp_R = lps[lp_num]['rprobe']
            Rs = np.append(Rs, tmp_R)

            # Average between times.
            lptimes = lps[lp_num]['time']
            idx_start = int(np.where(np.abs(lptimes-tmin) == np.abs(lptimes-tmin).min())[0][0])
            idx_end   = int(np.where(np.abs(lptimes-tmax) == np.abs(lptimes-tmax).min())[0][0])
            tmp_te = lps[lp_num]['temp'][idx_start:idx_end].mean()
            tmp_ne = lps[lp_num]['dens'][idx_start:idx_end].mean()
            tmp_V  = lps[lp_num]['pot'][idx_start:idx_end].mean()
            Tes = np.append(Tes, tmp_te)
            nes = np.append(nes, tmp_ne)
            Vs  = np.append(Vs,  tmp_V)

    # Remove unwanted LPs.
    if p_exclude:
        Rs  = np.delete(Rs, p_exclude)
        Tes = np.delete(Tes, p_exclude)
        nes = np.delete(nes, p_exclude)
        Vs  = np.delete(Vs, p_exclude)

    # Pressure as just nkT.
    ps = nes * Tes * 10**(-18)

    # Do an exponential fit.
    def exp_fit(x, a, b):
        return a * np.exp(-b * x)

    if V_adjust:
        Vs = Vs + np.abs(Vs.min()) + 0.01

    popt_V, pcov_V = curve_fit(exp_fit, Rs, Vs, maxfev=10000)
    popt_p, pcov_p = curve_fit(exp_fit, Rs, ps, maxfev=10000)
    popt_ne, pcov_ne = curve_fit(exp_fit, Rs, nes * 10**(-18), maxfev=10000)
    popt_te, pcov_te = curve_fit(exp_fit, Rs, Tes, maxfev=10000)

    # Get the average SP location.
    conn = MDSplus.Connection('localhost')
    ga_obj = gadata('rvsout', shot, connection=conn)
    ga_idx = np.logical_and(ga_obj.xdata > tmin, ga_obj.xdata < tmax)
    avg_rvsout = ga_obj.zdata[ga_idx].mean()

    # Average magnetic field.
    bt_obj = gadata('bt', shot, connection=conn)
    bt_idx = np.logical_and(bt_obj.xdata > tmin, bt_obj.xdata < tmax)
    avg_bt = bt_obj.zdata[bt_idx].mean()

    # Get derivative values at SP (assuming exponential so easy).
    sp_E     = popt_V[1] * exp_fit(avg_rvsout, *popt_V)
    sp_pgrad = popt_p[1] * exp_fit(avg_rvsout, *popt_p)
    sp_ne    = exp_fit(avg_rvsout, *popt_ne)
    sp_te    = exp_fit(avg_rvsout, *popt_te)

    # Print out info.
    print('Strike Point R:      {:.2f}'.format(avg_rvsout))
    print('Strike Point E:      {:.2e}'.format(sp_E))
    print('Strike Point grad p: {:.2e}'.format(sp_pgrad * 10**(24)))
    print('Strike Point ne:     {:.2e}'.format(sp_ne * 10**(24)))
    print('Strike Point Te:     {:.2f}'.format(sp_te))
    print('Average Bt:          {:.2f}'.format(avg_bt))

    fit_Rs = np.linspace(Rs.min(), Rs.max(), 100)
    fit_Vs = exp_fit(fit_Rs, *popt_V)
    fit_ps = exp_fit(fit_Rs, *popt_p)
    fit_nes = exp_fit(fit_Rs, *popt_ne)

    fig = pp.pplot(Rs, Vs)
    fig = pp.pplot(fit_Rs, fit_Vs, fig=fig, fmt='--', ylabel='V')

    fig = pp.pplot(Rs, ps*10**(24))
    fig = pp.pplot(fit_Rs, fit_ps*10**(24), fig=fig, fmt='--', ylabel='pe')

    fig = pp.pplot(Rs, nes*10**(6))
    fig = pp.pplot(fit_Rs, fit_nes*10**(24), fig=fig, fmt='--', ylabel='ne')


def use_ts(shot, tmin, tmax):
    """
    Calculate the PF stuff using the core TS system.
    """

    # Load in the TS data.
    ts = ThomsonClass(shot, 'core')
    ts.load_ts()
    ts.map_to_efit(np.linspace(tmin, tmax, 10))

    # Pressure as just nkT.
    Rs  = ts.avg_omp['RminRsep_omp']

    # Get the part only outside the separatrix.
    outside = np.where(Rs > 0.0)[0]
    Rs = Rs[outside]
    nes = ts.avg_omp['ne_omp'][outside] * 10**(-18)
    Tes = ts.avg_omp['Te_omp'][outside]
    ps = nes * Tes

    # Do an exponential fit.
    def exp_fit(x, a, b):
        return a * np.exp(-b * x)

    popt_p, pcov_p = curve_fit(exp_fit, Rs, ps, maxfev=10000)
    popt_ne, pcov_ne = curve_fit(exp_fit, Rs, nes, maxfev=10000)
    popt_te, pcov_te = curve_fit(exp_fit, Rs, Tes, maxfev=10000)

    fit_Rs = np.linspace(Rs.min(), Rs.max(), 100)
    fit_ps = exp_fit(fit_Rs, *popt_p)
    fit_nes = exp_fit(fit_Rs, *popt_ne)
    fit_tes = exp_fit(fit_Rs, *popt_te)

    fig = pp.pplot(Rs, ps)
    fig = pp.pplot(fit_Rs, fit_ps, fig=fig, fmt='--', ylabel='pe')

    fig = pp.pplot(Rs, nes)
    fig = pp.pplot(fit_Rs, fit_nes, fig=fig, fmt='--', ylabel='ne')

    fig = pp.pplot(Rs, Tes)
    fig = pp.pplot(fit_Rs, fit_tes, fig=fig, fmt='--', ylabel='Te')

    # dp/dr at the sep, multiply to make correct units.
    dpdr = -popt_p[1] * exp_fit(0, *popt_p) * 10**(18) * 1.609*10**(-19)
    nesep = exp_fit(0, *popt_ne) * 10**(18)
    print('dp/dr at separatrix (J/m3 or N/m2): {:.3e}'.format(dpdr))
    print('ne at separatrix (m-3):             {:.3e}'.format(nesep))
