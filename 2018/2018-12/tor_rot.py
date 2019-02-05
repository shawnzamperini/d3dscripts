import numpy as np
import pandas as pd
from gadata import gadata
import MDSplus as mds
from ThomsonClass import ThomsonClass  # For the load_gfile_mds function.
import pretty_plots as pp
from scipy.interpolate import Rbf, interp1d
import sys
from scipy.optimize import curve_fit


names = ['M18', 'V24', 'T41', 'M19', 'V09', 'T24', 'M20', 'V25', 'T10',
         'M21', 'V10', 'T42', 'M22', 'V26', 'T11', 'M23', 'V11', 'T43',
         'V27', 'M24', 'T12', 'M25', 'V12', 'T44', 'V28', 'M26', 'T13',
         'V13', 'M27', 'T45', 'V29', 'M28', 'T14', 'M29', 'V14', 'T46',
         'M30', 'V30', 'T15', 'M31', 'V15', 'T47', 'M32', 'V31', 'T16',
         'V16', 'T48', 'V32']

# Using 330RT radial values.
rads = [220.21, 219.64, 220.56, 220.93, 220.45, 221.3 , 221.72, 221.15,
        222.14, 222.46, 221.95, 222.77, 223.11, 222.72, 223.45, 223.79,
        223.5 , 224.13, 224.22, 224.6 , 224.9 , 225.27, 225.01, 225.62,
        225.69, 225.95, 226.34, 226.49, 226.71, 227.03, 227.21, 227.42,
        227.74, 228.08, 228.1 , 228.43, 228.82, 228.82, 229.2 , 229.44,
        229.55, 229.79, 230.21, 230.27, 230.62, 231.03, 231.13, 231.75]

names = np.array(names)
rads  = np.array(rads)

def run(shot, tmin=2500, tmax=4500, tmany=5, dist_outs=[0.5, 0.75], plot_it=True):
    """
    Map the (R, Z) of the 330RT CER chord to the OMP. Do an exponential fit
    to get the falloff. Get the toroidal rotation at the separatrix and
    probe tip (if possible).
    """

    times = np.linspace(tmin, tmax, tmany)
    rotc_df = pd.DataFrame(names).set_index(0)
    all_omps = np.zeros((tmany, len(names)))
    avg_rotcs = np.array([])

    omp_count = 0
    for time in times:

        # Load gfile.
        ts    = ThomsonClass(shot, 'core')
        gfile = ts.load_gfile_mds(shot, time)

        # The normal map to omp procedure. Put into cms.
        Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
        Rs = Rs * 100; Zs = Zs * 100
        Z_axis = gfile['ZmAxis'] * 100
        R_axis = gfile['RmAxis'] * 100
        Rs_trunc = Rs > R_axis

        # Only want the outboard half since thats where we're mapping R-Rsep OMP to.
        Zes_outboard = np.copy(gfile['lcfs'][:, 1][13:-17]) * 100
        Res_outboard = np.copy(gfile['lcfs'][:, 0][13:-17]) * 100

        # Interpolation functions of psin(R, Z) and R(psin, Z). Rs_trunc
        # helps with not interpolating the entire plasma, and just that
        # to the right of the magnetic axis, which is normally good enough.
        f_psiN = Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
        f_Romp = Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
        f_Rs = interp1d(Zes_outboard, Res_outboard, assume_sorted=False)

        rsep_omp = f_Rs(Z_axis)

        # Get the psin of each CER chord (which are all at Z=0).
        cer_psins = f_psiN(rads, np.full(len(rads), 0))

        # Then find the R value at the omp. Put in array holding all values.
        cer_omps = f_Romp(cer_psins, np.full(len(cer_psins), Z_axis))
        all_omps[omp_count] = cer_omps - rsep_omp

    avg_omps = np.mean(all_omps, axis=0)

    # Get the ROTC data.
    conn = mds.Connection('localhost')
    print('Loading chord: ', end='')
    for name in names:
        print(name, end=', ')
        tag = 'CERAROTC' + name
        try:
            gaobj = gadata(tag, shot, connection=conn, print_out=False)

            time = gaobj.xdata
            rotc = gaobj.zdata

            idx = np.logical_and(time>tmin, time<tmax)
            avg_rotc = np.mean(rotc[idx])
            avg_rotcs = np.append(avg_rotcs, avg_rotc)
        except:
            print('(!)', end='')

        #if name == 'T48':
        #    print('T48 ROTC: {:.2f}'.format(avg_rotc))

    # Find the first point inside the sep and outside the sep, do a line
    # between the two and then find what the value at r-rsep omp = 0 is. But
    # first, remove all the nans.
    no_nans   = ~np.isnan(avg_rotcs)
    avg_omps  = avg_omps[no_nans]
    avg_rotcs = avg_rotcs[no_nans]
    first_in  = np.where(avg_omps < 0)[0][-1]
    first_out = np.where(avg_omps > 0)[0][0]

    x1 = avg_omps[first_in]
    x2 = avg_omps[first_out]
    y1 = avg_rotcs[first_in]
    y2 = avg_rotcs[first_out]
    m = (y2 - y1)/(x2 - x1)

    # Point slope form, plugging in the point (0, y), where y = sep_rotc.
    # y - y1 = m * (x - x1)
    sep_rotc = m * (0 - x1) + y1

    print('\n\nROTC at separatrix: {:.2f}'.format(sep_rotc))
    print('(X1, Y1): ({:.2f}, {:.2f})'.format(x1, y1))
    print('(X2, Y2): ({:.2f}, {:.2f})'.format(x2, y2))
    print('Error: {:.3f}'.format(np.abs(y2-y1)/2))

    # Do the same but for a distance outside the separatrix.
    try:
        for dist_out in list(dist_outs):
            first_in  = np.where(avg_omps < dist_out)[0][-1]
            first_out = np.where(avg_omps > dist_out)[0][0]

            x1 = avg_omps[first_in]
            x2 = avg_omps[first_out]
            y1 = avg_rotcs[first_in]
            y2 = avg_rotcs[first_out]
            m = (y2 - y1)/(x2 - x1)

            print('(X1, Y1): ({:.2f}, {:.2f})'.format(x1, y1))
            print('(X2, Y2): ({:.2f}, {:.2f})'.format(x2, y2))

            dist_rotc = m * (dist_out - x1) + y1
            print('ROTC {:.2f} cm from separatrix: {:.2f}'.format(dist_out, dist_rotc))

    except:
        print('Could not get ROTC for {} cm from the separatrix.'.format(dist_out))

    #def exp_fit(x, a, b):
    #    return a * np.exp(-b * x)

    # Restrict the fit to just points outside the separatrix.
    #outside = 0 < avg_omps

    #clean_avg_omps  = avg_omps[outside][~np.isnan(avg_rotcs[outside])]
    #clean_avg_rotcs = avg_rotcs[outside][~np.isnan(avg_rotcs[outside])]
    #popt, pcov = curve_fit(exp_fit, clean_avg_omps, clean_avg_rotcs, maxfev=10000)
    #fit_omps  = np.linspace(0, avg_omps[outside].max(), 100)
    #fit_rotcs = exp_fit(fit_omps, *popt)

    #fig = pp.pplot(fit_omps, fit_rotcs, fmt='--')
    fig = pp.pplot([0], [sep_rotc], ms=12)
    fig = pp.pplot(avg_omps, avg_rotcs, fig=fig, fmt='.', xlabel='R-Rsep OMP (cm)',
                   ylabel='Toroidal Rotation (km/s)')

    # Plot the chord names over each data point.
    for i, name in enumerate(names[no_nans]):
        fig.axes[0].annotate(name, (avg_omps[i], avg_rotcs[i]))

    # Plot vertical lines at each dist_out.
    for dist in dist_outs:
        fig.axes[0].axvline(x=dist)

    return avg_omps, avg_rotcs
