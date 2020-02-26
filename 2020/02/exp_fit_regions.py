import numpy   as np
import tkinter as tk
import netCDF4
from scipy.optimize import curve_fit
from tkinter import filedialog
import matplotlib.pyplot as plt
import sys


# Colors and values for plotting.
red    = (214/255, 39/255,  40/255)
purple = (148/255, 103/255, 189/255)
fontsize = 16
lw = 2
ms = 12
plt.rcParams['font.family'] = 'DejaVu Sans'

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

def plot_with_fit(error_bands=True, ncpath=None, nc=None, log=True, region1_end=None,
                  fit=True, xlims=(0, 10), plot_a2=False, drop_tip=0, plot=1):
    """
    Plot centerline data from a 3DLIM run with optional exponential fits,
    error bands, and a comparision to A2 data. Currently just assumes the
    ITF side is divided into two regions, though extension to OTF would not
    be difficult.

    error_bands: Show error bands on plots, defined here as the standard
                  deviation of the 2D profile in the poloidal direction. It is
                  converted to a percent in the script.
    ncpath:      Path to the NetCDF file. If None, will open a file dialog.
    log:         Option to make the y-axis log or not.
    region1_end: Region 1 of the ITF side. Goes from 0-region1_end, and region 2
                  continues to the end of the probe. These are the regions where
                  two separate exponential fits are performed.
    fit:         Option to do exponential fit or not.
    xlims:       X-limits of the graph, passed as tuple, i.e. (low, high).
    plot_a2:     Option to plot A2 data as well for comparison.
    drop_itf_tip: An integer. Choose how many of the points on the tip of the
                   ITF you wish to drop. Sometimes there is an errant data point
                   that throws off the normalizations a bit, and it makes sense
                   to drop it.
    """

    if nc == None:
        if ncpath == None:
            # Select path to netcdf file.
            root = tk.Tk(); root.withdraw()
            initialdir = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/'
            ncpath = filedialog.askopenfilename(initialdir=initialdir,
                                                filetypes=(('NetCDF files', '*.nc'),),
                                                title='Select 3DLIM NetCDF file')
            print(ncpath)

        # Load in the NetCDF file.
        nc = netCDF4.Dataset(ncpath)

    # Load in the deposition array.
    dep_arr = np.array(nc.variables['NERODS3'][0] * -1)

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
    if error_bands:
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

    # Exponential fit to the different regions.
    if fit:

        def exp_fit(x, a, b):
            return a * np.exp(-b * x)

        # Define different regions for fit if applicable.
        if region1_end != None:
            region1_idx = np.where(itf_x <= region1_end)[0]
            region2_idx = np.where(itf_x >= region1_end)[0]
            popt1_itf, pcov1_itf = curve_fit(exp_fit, itf_x[region1_idx],
                                             itf_y[region1_idx], maxfev=5000)
            popt2_itf, pcov2_itf = curve_fit(exp_fit, itf_x[region2_idx],
                                             itf_y[region2_idx], maxfev=5000)
            region1_fitx = np.linspace(0, region1_end, 100)
            region2_fitx = np.linspace(region1_end, 15, 100)
            region1_fity_itf = exp_fit(region1_fitx, *popt1_itf)
            region2_fity_itf = exp_fit(region2_fitx, *popt2_itf)

            # Append together so we don't need to reuse this flag during plotting.
            fity_itf = np.append(region1_fity_itf, region2_fity_itf)
            fitx_itf = np.append(region1_fitx, region2_fitx)

            # Print out just so it's there.
            print("ITF")
            print("  Region 1: {:8.5f} * exp(-{:.3f} * x)".format(popt1_itf[0], 1/popt1_itf[1]))
            print("  Region 2: {:8.5f} * exp(-{:.3f} * x)".format(popt2_itf[0], 1/popt2_itf[1]))

        else:
            popt_itf, pcov_itf = curve_fit(exp_fit, itf_x, itf_y, maxfev=5000)
            fitx_itf = np.linspace(0, 15, 100)
            fity_itf = exp_fit(fitx_itf, *popt_itf)

        # Just assuming OTF has no regions, but easy enough to add in if we want.
        popt_otf, pcov_otf = curve_fit(exp_fit, otf_x, otf_y, maxfev=5000)
        fitx_otf = np.linspace(0, 15, 100)
        fity_otf = exp_fit(fitx_otf, *popt_otf)
        print("OTF")
        print("  Region 1: {:8.5f} * exp(-{:.3f} * x)".format(popt_otf[0], 1/popt_otf[1]))

    # Plotting commands.
    if 1 in plot:
        fig, ax = plt.subplots(figsize=(8,5))

        # Option to plot error bands or not.
        if error_bands:

            # Warnings are fine, just from nan's but they're handled correctly.
            old_settings = np.seterr(all='ignore')
            ax.fill_between(itf_x, itf_y-itf_y*itf_err, itf_y+itf_y*itf_err, color=red, alpha=0.35)
            ax.fill_between(otf_x, otf_y-otf_y*otf_err, otf_y+otf_y*otf_err, color=purple, alpha=0.35)
            np.seterr(**old_settings)

        # Plot exponential fits on top.
        if fit:
            ax.plot(fitx_itf, fity_itf, linestyle='--', color=red, lw=lw)
            ax.plot(fitx_otf, fity_otf, linestyle='--', color=purple, lw=lw)

            # Text labels for the fits.
            if region1_end != None:
                itf_str = r"$\lambda_W^1$ = {:.2f} cm".format(1/popt1_itf[1]) + "\n" + "$\lambda_W^2$ = {:.2f} cm".format(1/popt2_itf[1])
            else:
                itf_str = r"$\lambda_W$ = {:.2f} cm".format(1/popt_itf[1])
            otf_str = r"$\lambda_W$ = {:.2f} cm".format(1/popt_otf[1])
            ax.text(0.05, 0.15, itf_str, fontdict={'color':red, 'fontsize':fontsize}, transform=ax.transAxes)
            ax.text(0.05, 0.08, otf_str, fontdict={'color':purple, 'fontsize':fontsize}, transform=ax.transAxes)

        ax.plot(itf_x, itf_y, linestyle='-', color=red, label='ITF', lw=lw)
        ax.plot(otf_x, otf_y, linestyle='-', color=purple, label='OTF', lw=lw)

        # Include A2 data, normalized, for comparison.
        if plot_a2:
            max_a2_y = max(a2_itf_y.max(), a2_otf_y.max())
            a2_itf_y_norm = a2_itf_y / max_a2_y
            a2_otf_y_norm = a2_otf_y / max_a2_y
            ax.plot(a2_itf_x, a2_itf_y_norm, '.', ms=ms, color=red, markeredgecolor='k', markeredgewidth=1)
            ax.plot(a2_otf_x, a2_otf_y_norm, '.', ms=ms, color=purple, markeredgecolor='k', markeredgewidth=1)

        if log:
            ax.set_yscale('log')
        ax.set_ylim((5e-3, None))
        if xlims != None:
            ax.set_xlim(xlims)
        ax.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
        ax.set_ylabel("Deposition (normalized)", fontsize=fontsize)
        ax.legend(loc='upper right', fontsize=fontsize)
        fig.tight_layout()
        fig.show()

    if 2 in plot:

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5), sharex=True, sharey=True)

        # Option to plot error bands or not.
        if error_bands:
            old_settings = np.seterr(all='ignore')
            ax1.fill_between(itf_x, itf_y-itf_y*itf_err, itf_y+itf_y*itf_err, color=red, alpha=0.35)
            ax2.fill_between(otf_x, otf_y-otf_y*otf_err, otf_y+otf_y*otf_err, color=purple, alpha=0.35)
            np.seterr(**old_settings)

        # Plot exponential fits on top.
        if fit:
            ax1.plot(fitx_itf, fity_itf, linestyle='--', color=red, lw=lw)
            ax2.plot(fitx_otf, fity_otf, linestyle='--', color=purple, lw=lw)

            # Text labels for the fits.
            if region1_end != None:
                itf_str = r"$\lambda_W^1$ = {:.2f} cm".format(1/popt1_itf[1]) + "\n" + "$\lambda_W^2$ = {:.2f} cm".format(1/popt2_itf[1])
            else:
                itf_str = r"$\lambda_W$ = {:.2f} cm".format(1/popt_itf[1])
            otf_str = r"$\lambda_W$ = {:.2f} cm".format(1/popt_otf[1])
            ax1.text(0.6, 0.6, itf_str, fontdict={'color':red,    'fontsize':fontsize}, transform=ax1.transAxes)
            ax2.text(0.6, 0.6, otf_str, fontdict={'color':purple, 'fontsize':fontsize}, transform=ax2.transAxes)

        ax1.plot(itf_x, itf_y, linestyle='-', color=red, label='ITF', lw=lw)
        ax2.plot(otf_x, otf_y, linestyle='-', color=purple, label='OTF', lw=lw)

        # Include A2 data, normalized, for comparison.
        if plot_a2:
            max_a2_y = max(a2_itf_y.max(), a2_otf_y.max())
            a2_itf_y_norm = a2_itf_y / max_a2_y
            a2_otf_y_norm = a2_otf_y / max_a2_y
            ax1.plot(a2_itf_x, a2_itf_y_norm, '.', ms=ms, color=red, markeredgecolor='k', markeredgewidth=1)
            ax2.plot(a2_otf_x, a2_otf_y_norm, '.', ms=ms, color=purple, markeredgecolor='k', markeredgewidth=1)

        if log:
            ax1.set_yscale('log')
            ax2.set_yscale('log')
        ax1.set_ylim((5e-3, None))
        if xlims != None:
            ax1.set_xlim(xlims)
            ax2.set_xlim(xlims)
        ax1.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
        ax2.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
        ax1.set_ylabel("Deposition (normalized)", fontsize=fontsize)
        ax1.legend(loc='upper right', fontsize=fontsize)
        ax2.legend(loc='upper right', fontsize=fontsize)
        fig.tight_layout()
        fig.show()

    return nc

if __name__ == '__main__':

    if len(sys.argv) > 1:
        if sys.argv[1] == 'test':
            nc = plot_with_fit(ncpath='/mnt/c/Users/Shawn/Documents/' + \
                                 'd3d_work/3DLIM Runs/colprobe-z2-026c.nc',
                               region1_end=4, xlims=(0,10), plot_a2=True)
        else:
            print("Error: Supported options are 'test', ...")
    else:
        region1_end = float(input("Location of step (cm) in ITF (press enter if none): "))
        if region1_end == "":
            region1_end = None
            fit=False
        else:
            fit=True

        while True:
            ans = input("Plot A2 data (y/n)? ")
            if ans == 'y':
                plot_a2 = True
                break
            elif ans == 'n':
                plot_a2 = False
                break
            else:
                print("Error: Please choose 'y' or 'n'.")


        nc = plot_with_fit(fit=fit, region1_end=region1_end, plot_a2=plot_a2,
                           drop_tip=1, plot=[1, 2])

        xlims = None
        while True:
            ans = input("Tweak plot settings (y/n)? ")
            if ans == 'y':
                print("1: X-axis boundaries")
                setting = input("Choose setting: ")
                if setting == '1':
                    xlim_ans = input("Enter x-limits, separated by commas: ")
                    xlow  = float(xlim_ans.split(',')[0])
                    xhigh = float(xlim_ans.split(',')[1])
                    xlims = (xlow, xhigh)
                    other_ans = input("Change other settings (y/n)? ")
                    if other_ans == 'y':
                        continue
                    elif other_ans == 'n':
                        break
                    else:
                        print("Just going to assume you said yes.")
            elif ans == 'n':
                sys.exit()
            else:
                print("Please choose either 'y' or 'n'.")

        nc = plot_with_fit(fit=fit, region1_end=region1_end, plot_a2=plot_a2,
                           drop_tip=1, xlims=xlims, nc=nc, plot=[1, 2])
