import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4
from scipy.optimize import curve_fit


# Colors and values for plotting.
red    = (214/255, 39/255,  40/255)
purple = (148/255, 103/255, 189/255)
fontsize = 16
lw = 2
ms = 12
rsepx = True

# Data for A2 since it's common to compare results against it.
a2_itf_x = np.array([9.5, 9, 8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3,
                     2.5, 2, 1.5, 1, 0.5, 0])
a2_itf_y = np.array([0, 0, 0.001363767, 0.001363767, 0.0040913, 0, 0.002727533,
                     0.0040913, 0.002727533, 0.025911564, 0.006818833,
                     0.024547798, 0.087281059, 0.163651986, 0.265934478,
                     0.336850338, 0.377763335, 0.409129966, 0.475954527,
                     0.444587896])
a2_otf_x = np.array([9.1, 8.6, 8.1, 7.6, 7.1, 6.6, 6.1, 5.6, 5.1, 4.6, 4.1, 3.6,
                     3.1, 2.6, 2.1, 1.6, 1.1, 0.6])
a2_otf_y = np.array([0, 0.001363767, 0.002727533, 0.006818833, 0, 0.0040913,
                     0.0040913, 0.002727533, 0.008182599, 0.008182599,
                     0.013637666, 0.03545793, 0.040912997, 0.07500716,
                     0.085917293, 0.11455639, 0.143195488, 0.135012889])

# Normalize it.
max_a2_y = max(a2_itf_y.max(), a2_otf_y.max())
a2_itf_y_norm = a2_itf_y / max_a2_y
a2_otf_y_norm = a2_otf_y / max_a2_y

# Load in the NetCDF data.
nc1 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-051g.nc' # Best convective run.
nc2 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052a.nc' # Best diffusive run. Will be 052a I think.
#nc1 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-054a.nc' # Best convective run.
#nc2 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-053b.nc' # Best diffusive run. Will be 052a I think.

def get_centerline(ncpath, drop_tip):

    # Load in multiple NetCDF files, if they exist, and combine them into
    # one deposition array. The first files don't have a 1 in them, but the
    # laters can. Load in the NetCDF file.
    nc = netCDF4.Dataset(ncpath)

    # Load in the deposition array.
    dep_arr = np.array(nc.variables['NERODS3'][0] * -1)

    # Add on contributions from repeat runs.
    for i in range(1, 20):
        try:
            #print('Loading {}...'.format(i))
            ncpath_add = ncpath.split('.nc')[0] + str(i) + '.nc'
            nc = netCDF4.Dataset(ncpath_add)
            dep_arr = dep_arr + np.array(nc.variables['NERODS3'][0] * -1)
            print("Found additional run: {}".format(ncpath_add))
        except:
            pass

    # Location of each P bin, and its width.
    ps     = np.array(nc.variables['PS'][:].data)
    pwids  = np.array(nc.variables['PWIDS'][:].data)

    # Array of poloidal locations (i.e. the center of each P bin).
    pol_locs = ps - pwids/2.0

    # Distance cell centers along surface (i.e. the radial locations).
    rad_locs = np.array(nc.variables['ODOUTS'][:].data)

    # Get the centerline index (or closest to it).
    cline = np.abs(pol_locs).min()

    # Index the deposition array at the centerline for plotting.
    itf_x = rad_locs[np.where(rad_locs > 0.0)[0]] * 100
    itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
    otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1 * 100
    otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

    # Want to compute these and convert to a percent before we normalize the
    # data.
    #if error_bands:
    itf_err = dep_arr[:, np.where(rad_locs > 0.0)[0]].std(axis=0)
    otf_err = dep_arr[:, np.where(rad_locs < 0.0)[0]].std(axis=0)

    # Ignore division by zero errors.
    old_settings = np.seterr(divide='ignore', invalid='ignore')
    itf_err = itf_err / itf_y
    otf_err = otf_err / otf_y
    np.seterr(**old_settings)

    # Drop errant points at the tip if desired.
    itf_x = itf_x[drop_tip:]
    itf_y = itf_y[drop_tip:]
    itf_err = itf_err[drop_tip:]
    otf_x = otf_x[:-drop_tip]
    otf_y = otf_y[:-drop_tip]
    otf_err = otf_err[:-drop_tip]

    # Normalize the data so it can be compared to actual probe data.
    max_y = max(itf_y.max(), otf_y.max())
    itf_y = itf_y / max_y
    otf_y = otf_y / max_y

    return {'itf_x':itf_x, 'itf_y':itf_y, 'itf_err':itf_err, 'otf_x':otf_x,
           'otf_y':otf_y, 'otf_err':otf_err}

# If we want the xaxis to be R-Rsep OMP then apply the linear fit to
# go from location to R-Rsep. If not, just the set the parameters to
# values that won't affect anything.
if rsepx:
    rsepx1 = 1.0907
    rsepx2 = 6.9035
else:
    rsepx1 = 1.0
    rsepx2 = 0.0

# Exponential fit to the RBS data.
def exp_fit(x, a, b):
    return a * np.exp(-b * x)

# Fit the first 6 for region 1.
popt_itf1, pcov_itf1 = curve_fit(exp_fit, a2_itf_x[-6:], a2_itf_y_norm[-6:])
popt_itf2, pcov_itf2 = curve_fit(exp_fit, a2_itf_x[:-6], a2_itf_y_norm[:-6])
popt_otf1, pcov_otf1 = curve_fit(exp_fit, a2_otf_x[-6:], a2_otf_y_norm[-6:])
popt_otf2, pcov_otf2 = curve_fit(exp_fit, a2_otf_x[:-6], a2_otf_y_norm[:-6])

# Printout.
print("ITF")
print("  Lambda 1: {:.2f} cm".format(1/popt_itf1[1]))
print("  Lambda 2: {:.2f} cm".format(1/popt_itf2[1]))
print("OTF")
print("  Lambda 1: {:.2f} cm".format(1/popt_otf1[1]))
print("  Lambda 2: {:.2f} cm".format(1/popt_otf2[1]))

# Create arrays for the plot.
itf1_fit_x = np.linspace(a2_itf_x[-6:].min(), a2_itf_x[-6:].max(), 100)
itf2_fit_x = np.linspace(a2_itf_x[:-6].min(), a2_itf_x[:-6].max(), 100)
otf1_fit_x = np.linspace(a2_otf_x[-6:].min(), a2_otf_x[-6:].max(), 100)
otf2_fit_x = np.linspace(a2_otf_x[:-6].min(), a2_otf_x[:-6].max(), 100)
itf1_fit_y = exp_fit(itf1_fit_x, *popt_itf1)
itf2_fit_y = exp_fit(itf2_fit_x, *popt_itf2)
otf1_fit_y = exp_fit(otf1_fit_x, *popt_otf1)
otf2_fit_y = exp_fit(otf2_fit_x, *popt_otf2)

d1 = get_centerline(nc1, 5)
d2 = get_centerline(nc2, 5)

# Plotting commands.
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(10,5))

#ax1.plot(rsepx1*d1['itf_x']+rsepx2, d1['itf_y'], color=red,    lw=lw, label='Convective')
#ax1.plot(rsepx1*d2['itf_x']+rsepx2, d2['itf_y'], color=purple, lw=lw, label='Diffusive')

ax1.plot(rsepx1*d1['itf_x']+rsepx2, d1['itf_y'], color=red, linestyle=(0, (1, 1)), lw=lw, label='Convective')
ax1.plot(rsepx1*d2['itf_x']+rsepx2, d2['itf_y'], color=red, linestyle='-',  lw=lw, label='Diffusive')

ax2.plot(rsepx1*d1['otf_x']+rsepx2, d1['otf_y'], color=red,    lw=lw, label='Convective')
ax2.plot(rsepx1*d2['otf_x']+rsepx2, d2['otf_y'], color=purple, lw=lw, label='Diffusive')

ax1.plot(rsepx1*itf1_fit_x+rsepx2, itf1_fit_y, linestyle='--', color='k', lw=lw)
ax1.plot(rsepx1*itf2_fit_x+rsepx2, itf2_fit_y, linestyle='--', color='k', lw=lw)
ax2.plot(rsepx1*otf1_fit_x+rsepx2, otf1_fit_y,  linestyle='--', color='k', lw=lw)
ax2.plot(rsepx1*otf2_fit_x+rsepx2, otf2_fit_y,  linestyle='--', color='k', lw=lw)
ax1.plot(rsepx1*a2_itf_x+rsepx2, a2_itf_y_norm, '.', ms=ms, color=red,    markeredgecolor='k', markeredgewidth=1)
ax2.plot(rsepx1*a2_otf_x+rsepx2, a2_otf_y_norm, '.', ms=ms, color=purple, markeredgecolor='k', markeredgewidth=1)
ax1.set_yscale('log')
ax1.set_xlim([6.8, 14])
ax1.set_ylim([0.001, 2])
ax1.legend(fontsize=fontsize)
ax2.legend(fontsize=fontsize)
ax1.set_xlabel('R-Rsep OMP (cm)', fontsize=fontsize)
ax2.set_xlabel('R-Rsep OMP (cm)', fontsize=fontsize)
ax1.set_ylabel('Deposition (normalized)', fontsize=fontsize)
fig.tight_layout()
fig.show()
