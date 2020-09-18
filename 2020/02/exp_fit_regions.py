import numpy   as np
import tkinter as tk
import netCDF4
from scipy.optimize import curve_fit
from tkinter import filedialog
import matplotlib.pyplot as plt
import sys
from matplotlib.lines import Line2D


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
max_a2_y = max(a2_itf_y.max(), a2_otf_y.max())
a2_itf_y = a2_itf_y / max_a2_y
a2_otf_y = a2_otf_y / max_a2_y

# Files for Dperp scan.
#nc0_1 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-046c.nc'
#nc1   = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-046d.nc'
#nc5   = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-046e.nc'
#nc10  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-046a.nc'  # Should be same as 035.
#nc50  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-046f.nc'

# Files except runs without the fudge factor.
nc0_1 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052b.nc'
nc1   = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052c.nc'
nc5   = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052d.nc'
nc10  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052a.nc'  # Should be the best fit.
nc50  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052e.nc'

nc0_1 = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-060a.nc'
nc1   = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-060b.nc'
nc5   = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-060c.nc'
nc10  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052a.nc'  # Should be the best fit.
nc50  = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-052e.nc'

def get_centerline(ncpath, drop_tip):

    # Load in multiple NetCDF files, if they exist, and combine them into
    # one deposition array. The first files don't have a 1 in them, but the
    # laters can. Load in the NetCDF file.
    nc = netCDF4.Dataset(ncpath)

    # Load in the deposition array.
    dep_arr = np.array(nc.variables['NERODS3'][0] * -1)

    # Add on contributions from repeat runs.
    for i in range(1, 11):
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


def plot_with_fit(error_bands=True, ncpath=None, nc=None, log=True, region1_end=None,
                  fit=True, xlims=(0, 10), plot_a2=False, drop_tip=1, plot=[1,2],
                  drop_otf_fit=1, dperp_scan=False, rsepx=True):
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
    drop_tip:    An integer. Choose how many of the points on the tip of the
                   ITF you wish to drop. Sometimes there is an errant data point
                   that throws off the normalizations a bit, and it makes sense
                   to drop it.
    dperp_scan:  Show the results of scanning through Dperp.
    drop_otf_fit: Don't actually drop these points from the plot, just exclude
                   them from the fit.
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

    dict035 = get_centerline(ncpath, drop_tip)
    itf_x = dict035['itf_x']; itf_y = dict035['itf_y']; itf_err = dict035['itf_err']
    otf_x = dict035['otf_x']; otf_y = dict035['otf_y']; otf_err = dict035['otf_err']

    # Exponential fit to the different regions.
    if fit:

        def exp_fit(x, a, b):
            return a * np.exp(-b * x)

        # Define different regions for fit if applicable.
        if region1_end != None:
            region1_idx = np.where(a2_itf_x <= region1_end)[0]
            region2_idx = np.where(a2_itf_x >= region1_end)[0]
            popt1_itf, pcov1_itf = curve_fit(exp_fit, a2_itf_x[region1_idx],
                                             a2_itf_y[region1_idx], maxfev=5000)
            popt2_itf, pcov2_itf = curve_fit(exp_fit, a2_itf_x[region2_idx],
                                             a2_itf_y[region2_idx], maxfev=5000)
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
            popt_itf, pcov_itf = curve_fit(exp_fit, a2_itf_x, a2_itf_y, maxfev=5000)
            fitx_itf = np.linspace(0, 15, 100)
            fity_itf = exp_fit(fitx_itf, *popt_itf)

        # Just assuming OTF has no regions, but easy enough to add in if we want.
        popt_otf, pcov_otf = curve_fit(exp_fit, a2_otf_x[:-drop_otf_fit], a2_otf_y[:-drop_otf_fit], maxfev=5000)
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

        if dperp_scan:

            # Load centerline data for each.
            dict0_1 = get_centerline(nc0_1, drop_tip)
            dict1   = get_centerline(nc1, drop_tip)
            dict5   = get_centerline(nc5, drop_tip)
            dict50  = get_centerline(nc50, drop_tip)

            # If we want the xaxis to be R-Rsep OMP then apply the linear fit to
            # go from location to R-Rsep. If not, just the set the parameters to
            # values that won't affect anything.
            if rsepx:
                rsepx1 = 1.0907
                rsepx2 = 6.9035
            else:
                rsepx1 = 1.0
                rsepx2 = 0.0

            # Add to plot as light grey lines.
            ax1.plot(rsepx1*dict0_1['itf_x']+rsepx2, dict0_1['itf_y'], color='grey', alpha=0.5)
            ax1.plot(rsepx1*dict1['itf_x']+rsepx2,   dict1['itf_y'],   color='grey', alpha=0.5)
            ax1.plot(rsepx1*dict5['itf_x']+rsepx2,   dict5['itf_y'],   color='grey', alpha=0.5)
            ax1.plot(rsepx1*dict50['itf_x']+rsepx2,  dict50['itf_y'],  color='grey', alpha=0.5)
            ax2.plot(rsepx1*dict0_1['otf_x']+rsepx2, dict0_1['otf_y'], color='grey', alpha=0.5)
            ax2.plot(rsepx1*dict1['otf_x']+rsepx2,   dict1['otf_y'],   color='grey', alpha=0.5)
            ax2.plot(rsepx1*dict5['otf_x']+rsepx2,   dict5['otf_y'],   color='grey', alpha=0.5)
            ax2.plot(rsepx1*dict50['otf_x']+rsepx2,  dict50['otf_y'],  color='grey', alpha=0.5)

        # Option to plot error bands or not.
        if error_bands:
            old_settings = np.seterr(all='ignore')
            ax1.fill_between(rsepx1*itf_x+rsepx2, itf_y-itf_y*itf_err, itf_y+itf_y*itf_err, color=red, alpha=0.35)
            ax2.fill_between(rsepx1*otf_x+rsepx2, otf_y-otf_y*otf_err, otf_y+otf_y*otf_err, color=purple, alpha=0.35)
            np.seterr(**old_settings)

        # Plot exponential fits on top.
        if fit:
            ax1.plot(rsepx1*fitx_itf+rsepx2, fity_itf, linestyle='--', color="k", lw=lw)
            ax2.plot(rsepx1*fitx_otf+rsepx2, fity_otf, linestyle='--', color="k", lw=lw)

            # Text labels for the fits.
            if region1_end != None:
                itf_str = r"$\lambda_W^1$ = {:.2f} cm".format(1/popt1_itf[1]) + "\n" + "$\lambda_W^2$ = {:.2f} cm".format(1/popt2_itf[1])
            else:
                itf_str = r"$\lambda_W$ = {:.2f} cm".format(1/popt_itf[1])
            otf_str = r"$\lambda_W$ = {:.2f} cm".format(1/popt_otf[1])
            #ax1.text(0.6, 0.7, itf_str, fontdict={'color':red,    'fontsize':fontsize}, transform=ax1.transAxes)
            #ax2.text(0.6, 0.7, otf_str, fontdict={'color':purple, 'fontsize':fontsize}, transform=ax2.transAxes)

        ax1.plot(rsepx1*itf_x+rsepx2, itf_y, linestyle='-', color=red, label='ITF', lw=lw)
        ax2.plot(rsepx1*otf_x+rsepx2, otf_y, linestyle='-', color=purple, label='OTF', lw=lw)

        # Include A2 data, normalized, for comparison.
        if plot_a2:
            max_a2_y = max(a2_itf_y.max(), a2_otf_y.max())
            a2_itf_y_norm = a2_itf_y / max_a2_y
            a2_otf_y_norm = a2_otf_y / max_a2_y
            ax1.plot(rsepx1*a2_itf_x+rsepx2, a2_itf_y_norm, '.', ms=ms, color=red, markeredgecolor='k', markeredgewidth=1)
            ax2.plot(rsepx1*a2_otf_x+rsepx2, a2_otf_y_norm, '.', ms=ms, color=purple, markeredgecolor='k', markeredgewidth=1)

        if log:
            ax1.set_yscale('log')
            ax2.set_yscale('log')
        ax1.set_ylim((5e-3, None))
        if xlims != None:
            ax1.set_xlim(xlims)
            ax2.set_xlim(xlims)
        if rsepx:
            ax1.set_xlabel("R-Rsep OMP (cm)", fontsize=fontsize)
            ax2.set_xlabel("R-Rsep OMP (cm)", fontsize=fontsize)
        else:
            ax1.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
            ax2.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
        ax1.set_ylabel("Deposition (normalized)", fontsize=fontsize)

        # Custom legend.
        custom_lines1 = [Line2D([0], [0], color=red, lw=lw), Line2D([0], [0], color=red, ms=ms, marker='.', lw=0, mec='k')]
        custom_lines2 = [Line2D([0], [0], color=purple, lw=lw), Line2D([0], [0], color=purple, ms=ms, marker='.', lw=0, mec='k')]
        ax1.legend(custom_lines1, ['3DLIM', 'RBS'], fontsize=fontsize)
        ax2.legend(custom_lines2, ['3DLIM', 'RBS'], fontsize=fontsize)

        #ax1.legend(loc='upper right', fontsize=fontsize)
        #ax2.legend(loc='upper right', fontsize=fontsize)
        fig.tight_layout()
        fig.show()

    return nc

if __name__ == '__main__':

    if len(sys.argv) > 1:
        if sys.argv[1] == 'test':

            # With location as the x-axis.
            #nc = plot_with_fit(ncpath='/mnt/c/Users/Shawn/Documents/' + \
            #                     'd3d_work/3DLIM Runs/colprobe-z2-046a.nc',
            #                   region1_end=3, xlims=(0,6), plot_a2=True,
            #                   dperp_scan=True, drop_otf_fit=10, rsepx=False)

            # With R-Rsep OMP as the x-axis.
            #nc = plot_with_fit(ncpath='/mnt/c/Users/Shawn/Documents/' + \
            #                     'd3d_work/3DLIM Runs/colprobe-z2-046a.nc',
            #                   region1_end=3, plot_a2=True, xlims=(6.8, 13.8),
            #                   dperp_scan=True, drop_otf_fit=10, rsepx=True)

            # Updated to apply exponential fits to the RBS data.
            nc = plot_with_fit(ncpath='/mnt/c/Users/Shawn/Documents/' + \
                                 'd3d_work/3DLIM Runs/colprobe-z2-052a.nc',
                               region1_end=3, plot_a2=True, xlims=(6.8, 13.8),
                               dperp_scan=True, rsepx=True, plot=[2])
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
        drop_otf = 1
        dperp_scan = False
        while True:
            ans = input("Tweak plot settings (y/n)? ")
            if ans == 'y':
                print("1: X-axis boundaries")
                print("2: Ignore first X points in OTF fit")
                print("3: Show Dperp scan")
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
                elif setting == '2':
                    drop_otf = int(input("Number of data points from the tip to ignore: "))
                    other_ans = input("Change other settings (y/n)? ")
                    if other_ans == 'y':
                        continue
                    elif other_ans == 'n':
                        break
                    else:
                        print("Just going to assume you said yes.")
                elif setting == '3':
                    dperp_scan=True
                    other_ans = input("Change other settings (y/n)? ")
                    if other_ans == 'y':
                        continue
                    elif other_ans == 'n':
                        break
            elif ans == 'n':
                sys.exit()
            else:
                print("Please choose either 'y' or 'n'.")

        nc = plot_with_fit(fit=fit, region1_end=region1_end, plot_a2=plot_a2,
                           drop_tip=1, xlims=xlims, nc=nc, plot=[1, 2],
                           drop_otf=drop_otf, dperp_scan=dperp_scan)
