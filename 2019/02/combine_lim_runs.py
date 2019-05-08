import netCDF4
import numpy as np
import pretty_plots as pp
from scipy.optimize import curve_fit

lim53_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z1-015c1.nc'
lim69_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z1-015c2.nc'

lim53 = netCDF4.Dataset(lim53_path)
lim69 = netCDF4.Dataset(lim69_path)

def centerline(lim1, lim2, lim2_mult):

    #The deposition array.
    dep_arr1 = np.array(lim1.variables['NERODS3'][0] * -1)
    dep_arr2 = np.array(lim2.variables['NERODS3'][0] * -1)

    # Location of each P bin, and its width. Currently they all have the same width,
    # but it may end up such that there are custom widths so we leave it like this.
    ps1     = np.array(lim1.variables['PS'][:].data)
    ps2     = np.array(lim2.variables['PS'][:].data)
    pwids1  = np.array(lim1.variables['PWIDS'][:].data)
    pwids2  = np.array(lim2.variables['PWIDS'][:].data)

    # Array of poloidal locations (i.e. the center of each P bin).
    pol_locs1 = ps1 - pwids1/2.0
    pol_locs2 = ps2 - pwids2/2.0

    # Drop last row since it's garbage.
    dep_arr1 = dep_arr1[:-1, :]
    dep_arr2 = dep_arr2[:-1, :]
    pol_locs1 = pol_locs1[:-1]
    pol_locs2 = pol_locs2[:-1]

    # Distance cell centers along surface (i.e. the radial locations).
    rad_locs1 = np.array(lim1.variables['ODOUTS'][:].data)
    rad_locs2 = np.array(lim2.variables['ODOUTS'][:].data)

    # Get the centerline index (or closest to it).
    cline1 = np.abs(pol_locs1).min()
    cline2 = np.abs(pol_locs2).min()

    # Index the deposition array at the centerline for plotting.
    itf_x1 = rad_locs1[np.where(rad_locs1 > 0.0)[0]]
    itf_y1 = dep_arr1[np.where(pol_locs1 == cline1)[0], np.where(rad_locs1 > 0.0)[0]]
    otf_x1 = rad_locs1[np.where(rad_locs1 < 0.0)[0]] * -1
    otf_y1 = dep_arr1[np.where(pol_locs1 == cline1)[0], np.where(rad_locs1 < 0.0)[0]]
    itf_x2 = rad_locs2[np.where(rad_locs2 > 0.0)[0]]
    itf_y2 = dep_arr2[np.where(pol_locs2 == cline2)[0], np.where(rad_locs2 > 0.0)[0]]
    otf_x2 = rad_locs2[np.where(rad_locs2 < 0.0)[0]] * -1
    otf_y2 = dep_arr2[np.where(pol_locs2 == cline2)[0], np.where(rad_locs2 < 0.0)[0]]

    itf_tot = itf_y1 + itf_y2 * lim2_mult
    otf_tot = otf_y1 + otf_y2 * lim2_mult

    fig = pp.pplot(itf_x1, itf_tot, fmt='-', label='ITF')
    fig = pp.pplot(otf_x1, otf_tot, fmt='-', fig=fig, color=8)

    def exp_fit(x, a, b):
        return a * np.exp(-b * x)

    popt_itf, pcov_itf = curve_fit(exp_fit, itf_x1, itf_tot, maxfev=5000)
    popt_otf, pcov_otf = curve_fit(exp_fit, otf_x1, otf_tot, maxfev=5000)

    itf_x_fit = np.linspace(0.01, 0.14, 100)
    itf_y_fit = exp_fit(itf_x_fit, *popt_itf)
    otf_x_fit = np.linspace(0.01, 0.14, 100)
    otf_y_fit = exp_fit(otf_x_fit, *popt_otf)

    fig = pp.pplot(itf_x_fit, itf_y_fit, fmt='--', fig=fig)
    fig = pp.pplot(otf_x_fit, otf_y_fit, fmt='--', fig=fig, color=8,
                   xlabel='Distance along probe (m)',
                   ylabel='W Deposition (arbitrary units)',
                   label='OTF')

    print("\na * exp(-b * x):")
    print("  ITF: a={:5.2f}  b={:5.2f}  1/b={:5.2f} cm".format(popt_itf[0], popt_itf[1], 1/popt_itf[1] * 100))
    print("  OTF: a={:5.2f}  b={:5.2f}  1/b={:5.2f} cm".format(popt_otf[0], popt_otf[1], 1/popt_otf[1] * 100))
    print("Max ITF/OTF Ratio: {:.2f}".format(itf_tot.max() / otf_tot.max()))
    print("Total ITF/OTF Ratio: {:.2f}".format(itf_tot.sum() / otf_tot.sum()))

def deposition_contour(lim1, lim2, side, lim2_mult=1.0, rad_cutoff=0.05):

    #The deposition array.
    dep_arr1 = np.array(lim1.variables['NERODS3'][0] * -1)
    dep_arr2 = np.array(lim2.variables['NERODS3'][0] * -1)

    # Location of each P bin, and its width. Currently they all have the same width,
    # but it may end up such that there are custom widths so we leave it like this.
    ps1     = np.array(lim1.variables['PS'][:].data)
    ps2     = np.array(lim2.variables['PS'][:].data)
    pwids1  = np.array(lim1.variables['PWIDS'][:].data)
    pwids2  = np.array(lim2.variables['PWIDS'][:].data)

    # Array of poloidal locations (i.e. the center of each P bin).
    pol_locs1 = ps1 - pwids1/2.0
    pol_locs2 = ps2 - pwids2/2.0

    # Drop last row since it's garbage.
    dep_arr1 = dep_arr1[:-1, :]
    dep_arr2 = dep_arr2[:-1, :]
    pol_locs1 = pol_locs1[:-1]
    pol_locs2 = pol_locs2[:-1]

    # Distance cell centers along surface (i.e. the radial locations).
    rad_locs1 = np.array(lim1.variables['ODOUTS'][:].data)
    rad_locs2 = np.array(lim2.variables['ODOUTS'][:].data)

    # Remove data beyond rad_cutoff.
    idx = np.where(np.abs(rad_locs1)<rad_cutoff)[0]
    rad_locs1 = rad_locs1[idx]
    rad_locs2 = rad_locs2[idx]
    dep_arr1 = dep_arr1[:, idx]
    dep_arr2 = dep_arr2[:, idx]

    # Get only positive values of rad_locs for ITF...
    idx = np.where(rad_locs1 > 0.0)[0]
    X_itf1, Y_itf1 = np.meshgrid(rad_locs1[idx], pol_locs1)
    X_itf2, Y_itf2 = np.meshgrid(rad_locs2[idx], pol_locs2)
    Z_itf1 = dep_arr1[:, idx]
    Z_itf2 = dep_arr2[:, idx]

    # ... negative for OTF.
    idx = np.where(rad_locs1 < 0.0)[0]
    X_otf1, Y_otf1 = np.meshgrid(np.abs(rad_locs1[idx][::-1]), pol_locs1)
    X_otf2, Y_otf2 = np.meshgrid(np.abs(rad_locs2[idx][::-1]), pol_locs2)
    Z_otf1 = dep_arr1[:, idx][:, ::-1]
    Z_otf2 = dep_arr2[:, idx][:, ::-1]

    # Make the levels for the contour plot out of whichever side has the max deposition.
    #if Z_itf.max() > Z_otf.max():
    #    levels = np.linspace(0, Z_itf.max(), 15)
    #else:
    #    levels = np.linspace(0, Z_otf.max(), 15)

    # Plotting commands.
    if side == 'ITF':
        X = X_itf1; Y = Y_itf1; Z = Z_itf1 + Z_itf2 * lim2_mult
    else:
        X = X_otf1; Y = Y_otf1; Z = Z_otf1 + Z_itf2 * lim2_mult

    fig = pp.pcontourf(X, Y, Z, cbarlabel=side, yrange=[-0.015, 0.015])



# Can this multiplier be used as a guage to look for trend in: higher mult (=higher sheet)
# would mean something more in the quantified peaking on the edges? Maybe change
# this variable to study ITF/OTF effect?
lim2_mult = 1.0
centerline(lim53, lim69, lim2_mult)
deposition_contour(lim53, lim69, 'ITF', lim2_mult)
deposition_contour(lim53, lim69, 'OTF', lim2_mult)
