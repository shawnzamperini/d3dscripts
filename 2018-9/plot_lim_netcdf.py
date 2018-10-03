import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pretty_plots as pp
from scipy.optimize import curve_fit


test = '62'
probe = 'b8'
filename = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-'+probe+'_test'+test+'.nc'
net = netCDF4.Dataset(filename)
print('Test: ' + test)

# 2D array of deposition data.
dep_arr = np.array(net.variables['NERODS3'][0] * -1)

# Location of each P bin, and its width. Currently they all have the same width,
# but it may end up such that there are custom widths so we leave it like this.
ps     = np.array(net.variables['PS'][:].data)
pwids  = np.array(net.variables['PWIDS'][:].data)

# Array of poloidal locations (i.e. the center of each P bin).
pol_locs = ps - pwids/2.0

# Distance cell centers along surface (i.e. the radial locations).
rad_locs = np.array(net.variables['ODOUTS'][:].data)

# Drop last row since it's garbage.
dep_arr = dep_arr[:-1, :]
pol_locs = pol_locs[:-1]

itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
itf_y = dep_arr[np.where(pol_locs == 0)[0], np.where(rad_locs > 0.0)[0]]
otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
otf_y = dep_arr[np.where(pol_locs == 0)[0], np.where(rad_locs < 0.0)[0]]

def plot_3d():
    Y, X = np.meshgrid(rad_locs, pol_locs)
    font = {'fontsize':18}
    fig = plt.figure(figsize=((15,8)))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, dep_arr, cmap='Reds')
    ax.set_xlabel('\nPoloidal (m)', font)
    ax.set_ylabel('\nDistance along probe (m)', font)
    ax.set_zlabel('\nDeposition counts', font)
    #fig.tight_layout()
    fig.show()

def plot_centerline():
    # Do centerline plots.
    fig = pp.pplot(itf_x, itf_y, fmt='-', label='ITF')
    fig = pp.pplot(otf_x, otf_y, fmt='-', fig=fig, color=8, label='OTF',
                   xlabel='Distance along probe (m)',
                   ylabel='Deposition (arbitrary units)')

def lambdas():
    def exp_fit(x, a, b):
        return a * np.exp(-x*b)

    itf_popt, itf_pcov = curve_fit(exp_fit, itf_x, itf_y)
    otf_popt, otf_pcov = curve_fit(exp_fit, otf_x, otf_y)
    itf_x_fit = np.linspace(itf_x.min(), itf_x.max(), 100)
    otf_x_fit = np.linspace(otf_x.min(), otf_x.max(), 100)
    itf_y_fit = exp_fit(itf_x_fit, *itf_popt)
    otf_y_fit = exp_fit(otf_x_fit, *otf_popt)

    print('ITF Lambda = {:.4f}'.format(1/itf_popt[1]))
    print('OTF Lambda = {:.4f}'.format(1/otf_popt[1]))

    fig = None
    fig = pp.pplot(itf_x_fit, itf_y_fit, '--', label='ITF Fit', color=6, fig=fig)
    fig = pp.pplot(otf_x_fit, otf_y_fit, '--', label='OTF Fit', color=8, fig=fig)
    fig = pp.pplot(itf_x, itf_y, '.', label='ITF', fig =fig)
    fig = pp.pplot(otf_x, otf_y, '.', label='OTF', color=8, fig=fig, xlabel='Distance along probe (m)', ylabel='Deposition Counts')

def itf_otf_content():
    maxw_ratio = max(itf_y) / max(otf_y)
    totw_ratio = sum(itf_y) / sum(otf_y)

    print('ITF/OTF Max W: {:.4f}'.format(maxw_ratio))
    print('ITF/OTF Tot W: {:.4f}'.format(totw_ratio))

def avg_pol_prof():
    print("Poloidal locations:")
    for pol in pol_locs:
        print(pol)
    print("Average values:")
    for val in np.mean(dep_arr, axis=1):
        print(val)

#plot_centerline()
#lambdas()
#itf_otf_content()
#plot_3d()
avg_pol_prof()
