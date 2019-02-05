import gadata
import numpy             as np
import get_ts            as ts
import numpy             as np
import matplotlib.pyplot as plt
import MDSplus           as mds
from scipy.optimize      import curve_fit


def fit_ne(shot, lowlim=None, uplim=None, tmin=2000, tmax=4000):
    ts_dict = ts.run_script(shot, 'core', tmin=tmin, tmax=tmax)

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

    def exp_fit(x, a, b):
        return a * np.exp(-x*b)

    popt_te,   pcov_te     = curve_fit(exp_fit, x_ts, y_te)
    popt_ne,   pcov_ne     = curve_fit(exp_fit, x_ts, y_ne)
    x_fit            = np.linspace(0, 14, 100)
    y_fit_te         = exp_fit(x_fit, *popt_te)
    y_fit_ne         = exp_fit(x_fit, *popt_ne)

    print("ne fall off length:       {:.3f}".format(1/popt_ne[1]))
    errs = np.sqrt(np.diag(pcov_ne))
    print("ne fall off length error: {:.4f}".format((errs[1]/popt_ne[1]) * (1/popt_ne[1])))
    print("Te fall off length:       {:.3f}".format(1/popt_te[1]))
    errs = np.sqrt(np.diag(pcov_te))
    print("Te fall off length error: {:.4f}".format((errs[1]/popt_te[1]) * (1/popt_te[1])))

    font = {'fontsize':24, 'weight':'bold'}
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


def get_avg_params(conn, shot, lower_time=2000, upper_time=4000):
    ne_dict = {}; ne_err_dict = {}
    a_dict  = {}; a_err_dict  = {}
    ip_dict = {}; ip_err_dict = {}
    bt_dict = {}; bt_err_dict = {}
    q95_dict = {}; q95_err_dict = {}
    pinj_dict = {}; pinj_err_dict = {}

    print("Getting values for shot: {}".format(shot))

    ne_obj  = gadata.gadata('\DENSITY', shot, connection=conn)
    times   = ne_obj.xdata
    nes     = ne_obj.zdata
    idx     = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_nes = np.mean(nes[idx]) * 100**(3) * 10**(-19)  # convert to m-3 then to 1E19 m-3
    std_nes = np.std(nes[idx]) * 100**(3) * 10**(-19)

    a_obj     = gadata.gadata('\AMINOR', shot, connection=conn)
    times     = a_obj.xdata
    amins     = a_obj.zdata
    idx       = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_amins = np.mean(amins[idx])
    std_amins = np.std(amins[idx])

    ip_obj  = gadata.gadata('\IP', shot, connection=conn)
    times   = ip_obj.xdata.data()
    ips     = ip_obj.zdata.data()
    idx     = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_ips = np.mean(ips[idx]) / 10**6   # convert to MA
    std_ips = np.std(ips[idx]) / 10**6

    q95_obj  = gadata.gadata('\Q95', shot, connection=conn)
    times = q95_obj.xdata
    q95s  = q95_obj.zdata
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_q95s = np.mean(q95s[idx])
    std_q95s = np.std(q95s[idx])

    bt_obj  = gadata.gadata('\BT', shot, connection=conn)
    times = bt_obj.xdata.data()
    bts  = bt_obj.zdata.data()
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_bts = np.mean(bts[idx])
    std_bts = np.std(bts[idx])

    pinj_obj  = gadata.gadata('\PINJ', shot, connection=conn)
    times = pinj_obj.xdata
    pinjs  = pinj_obj.zdata
    idx    = np.where(np.logical_and(times>=lower_time, times<=upper_time))
    avg_pinjs = np.mean(pinjs[idx])
    std_pinjs = np.std(pinjs[idx])

    print("  Line Avg. Density (1E19 m-3): {:.3f} +- {:.4f}".format(avg_nes, std_nes))
    print("  Minor Radius (m): {:.3f} +- {:.4f}".format(avg_amins, std_amins))
    print("  Plasma Current (MA): {:.3f} +- {:.4f}".format(avg_ips, std_ips))
    print("  Q95: {:.3f} +- {:.4f}".format(avg_q95s, std_q95s))
    print("  Toroidal Mag. (T): {:.3f} +- {:.4f}".format(avg_bts, std_bts))
    print("  Power Injected (MW): {:.3f} +- {:.4f}".format(avg_pinjs, std_pinjs))
