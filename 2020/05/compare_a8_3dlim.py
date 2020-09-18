import lim_plots
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


# Some constants.
itf_cal = 0.5E-06  # LAMS calibration to go from counts to areal density.
otf_cal = 0.5E-06
tip_itf_ignore = 3  # Ignore the first this many points on 3DLIM.
tip_otf_ignore = 3
tip_itf_lams_ignore = 18  # Ignore the first this many point on LAMS.
tip_otf_lams_ignore = 18
log_plot       = True
itfotf_plot    = True
plot_opt       = 2
lim_band       = False
lim_band_width = 5
lam_band       = False
smooth         = True

# Path to the 3DLIM data.
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-045h.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-045e.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-001v.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-001w.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-002b.nc'
ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-004u.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-005a.nc'

# My colors.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (0, 0, 0)]
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# Path to A8 RBS data.
a8_rbs = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
         'Collector Probe Excel Sheets/A8.xlsx'

# Path to A8 LAMS data.
au8_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AU08_Map_Analysis.xlsx'
ad8_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AD08_Map_Analysis.xlsx'

# Load 3DLIM data.
lp = lim_plots.LimPlots(ncpath)

# Load into Dataframes.
rbs = pd.read_excel(a8_rbs, nrows=20)
au8_lams = pd.read_excel(au8_file)
ad8_lams = pd.read_excel(ad8_file)

# Pivot the data so the index is the poloidal data and the columns the radial.
au8_lams = au8_lams.pivot(columns='z Location [mm]', index='Axial Location [mm]')
ad8_lams = ad8_lams.pivot(columns='z Location [mm]', index='Axial Location [mm]')

# Pull out the data we will plot. RBS first.
rbs_itf_x     = rbs['Distance from Tip D (cm)'].values
rbs_itf_y     = rbs['W Areal Density D (1e15 W/cm2)'].values
rbs_itf_y_err = rbs['W Areal Density Error D (1e15 W/cm2)'].values
rbs_otf_x     = rbs['Distance from Tip U (cm)'].values
rbs_otf_y     = rbs['W Areal Density U (1e15 W/cm2)'].values
rbs_otf_y_err = rbs['W Areal Density Error U (1e15 W/cm2)'].values

# Then the LAMS values of the centerline data.
#lams_itf_x     = ad8_lams.index.values / 10  # mm to cm
#lams_itf_y     = ad8_lams.mean(axis=1).values
#lams_otf_x     = au8_lams.index.values / 10  # mm to cm
#lams_otf_y     = au8_lams.mean(axis=1).values
#lams_itf_y_err = ad8_lams.std(axis=1).values
#lams_otf_y_err = au8_lams.std(axis=1).values

# Get the centerline data.
lams_itf_x = ad8_lams["Total W"][2.5].index.values / 10  # mm to cm
lams_itf_y = ad8_lams["Total W"][2.5].values
lams_itf_y_err = ad8_lams["Total W error"][2.5].values
lams_otf_x = au8_lams["Total W"][2.5].index.values / 10  # mm to cm
lams_otf_y = au8_lams["Total W"][2.5].values
lams_otf_y_err = au8_lams["Total W error"][2.5].values

# LAMS data of the data + err and - err so we can apply all the normalization
# and data cleanup to it as well.
lams_itf_y_errp = lams_itf_y + lams_itf_y_err
lams_itf_y_errm = lams_itf_y - lams_itf_y_err
lams_otf_y_errp = lams_otf_y + lams_otf_y_err
lams_otf_y_errm = lams_otf_y - lams_otf_y_err

# Apply calibrations. Won't matter in the end if the calibrations are the same,
# since it all gets normalized, but it will if they're different.
lams_itf_y     *= itf_cal
lams_itf_y_err *= itf_cal
lams_otf_y     *= otf_cal
lams_otf_y_err *= otf_cal

# And then the 3DLIM data.
lim_dict = lp.centerline(show_plot=False)
lim_itf_x = lim_dict['itf_x'] * 100  # m to cm
lim_itf_y = lim_dict['itf_y']
lim_otf_x = lim_dict['otf_x'] * 100  # m to cm
lim_otf_y = lim_dict['otf_y']

# Shift the LAMS data so the minimum value is zero.
min_lams = min(lams_itf_y.min(), lams_otf_y.min())
lams_itf_y = lams_itf_y - min_lams
lams_otf_y = lams_otf_y - min_lams

# Reverse the OTF data so it goes from tip to end.
lim_otf_x = lim_otf_x[::-1]
lim_otf_y = lim_otf_y[::-1]

# Ignore number of data points.
lim_itf_x = lim_itf_x[tip_itf_ignore:]
lim_itf_y = lim_itf_y[tip_itf_ignore:]
lim_otf_x = lim_otf_x[tip_otf_ignore:]
lim_otf_y = lim_otf_y[tip_otf_ignore:]
lams_itf_x = lams_itf_x[tip_itf_lams_ignore:]
lams_itf_y = lams_itf_y[tip_itf_lams_ignore:]
lams_otf_x = lams_otf_x[tip_otf_lams_ignore:]
lams_otf_y = lams_otf_y[tip_otf_lams_ignore:]
lams_itf_y_errp = lams_itf_y_errp[tip_itf_lams_ignore:]
lams_itf_y_errm = lams_itf_y_errm[tip_itf_lams_ignore:]
lams_otf_y_errp = lams_otf_y_errp[tip_otf_lams_ignore:]
lams_otf_y_errm = lams_otf_y_errm[tip_otf_lams_ignore:]

# Interpolate over some spurious data points.
def replace_linear(arr_x, arr_y, start_idx, end_idx):
    x1 = arr_x[start_idx]
    x2 = arr_x[end_idx]
    y1 = arr_y[start_idx]
    y2 = arr_y[end_idx]
    m = (y2 - y1) / (x2 - x1)
    arr_y[start_idx:end_idx+1] = m * (arr_x[start_idx:end_idx+1] - x1) + y1
    return arr_y

def replace_exp(arr_x, arr_y, start, end):
    def exp_fit(x, a, b):
        return a * np.exp(-b * x)

    # Find indicies of the regions.
    start_idx = np.where(np.abs(arr_x-start) == np.min(np.abs(arr_x-start)))[0][0]
    end_idx = np.where(np.abs(arr_x-end) == np.min(np.abs(arr_x-end)))[0][0]
    #print(start_idx)
    #print(end_idx)
    #print("Old: {}".format(arr_y[start_idx:end_idx]))

    # Make the fit out of the bounding two data points.
    for_fit_x = np.append(arr_x[start_idx-2:start_idx], arr_x[end_idx:end_idx+2])
    for_fit_y = np.append(arr_y[start_idx-2:start_idx], arr_y[end_idx:end_idx+2])
    popt, pcov = curve_fit(exp_fit, for_fit_x, for_fit_y)
    #print(popt)

    # Replace with fit values.
    arr_y[start_idx:end_idx] = exp_fit(arr_x[start_idx:end_idx], *popt)
    #print("New: {}".format(arr_y[start_idx:end_idx]))
    return arr_y

lams_itf_y[21] = (lams_itf_y[22] + lams_itf_y[20]) / 2
lams_itf_y[281] = (lams_itf_y[280] + lams_itf_y[282]) / 2

"""
# Fix some errant data.
lams_itf_y[204] = (lams_itf_y[203] + lams_itf_y[205]) / 2
lams_itf_y_errp[204] = (lams_itf_y_errp[203] + lams_itf_y_errp[205]) / 2
lams_itf_y_errm[204] = (lams_itf_y_errm[203] + lams_itf_y_errm[205]) / 2
lams_otf_y = replace_exp(lams_otf_x, lams_otf_y, 0.61, 1.00)
lams_otf_y = replace_exp(lams_otf_x, lams_otf_y, 2.13, 2.60)
lams_otf_y = replace_exp(lams_otf_x, lams_otf_y, 3.10, 3.67)
lams_otf_y = replace_exp(lams_otf_x, lams_otf_y, 5.12, 5.21)
lams_otf_y_errp = replace_exp(lams_otf_x, lams_otf_y_errp, 0.61, 1.00)
lams_otf_y_errp = replace_exp(lams_otf_x, lams_otf_y_errp, 2.13, 2.60)
lams_otf_y_errp = replace_exp(lams_otf_x, lams_otf_y_errp, 3.10, 3.67)
lams_otf_y_errp = replace_exp(lams_otf_x, lams_otf_y_errp, 5.12, 5.21)
lams_otf_y_errm = replace_exp(lams_otf_x, lams_otf_y_errm, 0.61, 1.00)
lams_otf_y_errm = replace_exp(lams_otf_x, lams_otf_y_errm, 2.13, 2.60)
lams_otf_y_errm = replace_exp(lams_otf_x, lams_otf_y_errm, 3.10, 3.67)
lams_otf_y_errm = replace_exp(lams_otf_x, lams_otf_y_errm, 5.12, 5.21)
"""

# Normalize all the data.
max_rbs = max(rbs_itf_y.max(), rbs_otf_y.max())
rbs_itf_y = rbs_itf_y / max_rbs
rbs_otf_y = rbs_otf_y / max_rbs
max_lams = max(lams_itf_y.max(), lams_otf_y.max())
lams_itf_y = lams_itf_y / max_lams
lams_otf_y = lams_otf_y / max_lams
lams_itf_y_errp = lams_itf_y_errp / max_lams
lams_itf_y_errm = lams_itf_y_errm / max_lams
lams_otf_y_errp = lams_otf_y_errp / max_lams
lams_otf_y_errm = lams_otf_y_errm / max_lams
max_lim = max(lim_itf_y.max(), lim_otf_y.max())
lim_itf_y = lim_itf_y / max_lim
lim_otf_y = lim_otf_y / max_lim

# Smooth out the data.
if smooth:
    window = 21
    lim_otf_y = savgol_filter(lim_otf_y, window, 3)
    lim_itf_y = savgol_filter(lim_itf_y, window, 3)

# ITF/OTF plot: Need to interpolate onto common x values to calculate ratios.
num_points = len(lams_itf_x)
lams_com_x = np.linspace(max(lams_itf_x.min(), lams_otf_x.min()), min(lams_itf_x.max(), lams_otf_x.max()), num_points)
lim_com_x  = np.linspace(max(lim_itf_x.min(),  lim_otf_x.min()),  min(lim_itf_x.max(),  lim_otf_x.max()),  num_points)
f_lams_itf = interp1d(lams_itf_x, lams_itf_y)
f_lams_otf = interp1d(lams_otf_x, lams_otf_y)
f_lim_itf  = interp1d(lim_itf_x,  lim_itf_y)
f_lim_otf  = interp1d(lim_otf_x,  lim_otf_y)
lams_itfotf = f_lams_itf(lams_com_x) / f_lams_otf(lams_com_x)
lim_itfotf  = f_lim_itf(lim_com_x)   / f_lim_otf(lim_com_x)

# Calculate arrays for 3DLIM band plot.
if lim_band:

    # First need to make sure the array we're using is divisible by our number
    # points to average.
    num_keep = lim_band_width*int(len(lim_itf_x) / lim_band_width)
    lim_itf_x_band = lim_itf_x[:num_keep]
    lim_itf_x_band = np.mean(lim_itf_x_band.reshape(-1, lim_band_width), axis=1)
    lim_itf_y_band = lim_itf_y[:num_keep]
    lim_itf_y_band_avg = np.mean(lim_itf_y_band.reshape(-1, lim_band_width), axis=1)
    lim_itf_y_band_std = np.std(lim_itf_y_band.reshape(-1, lim_band_width), axis=1)

    # OTF side.
    num_keep = lim_band_width*int(len(lim_otf_x) / lim_band_width)
    lim_otf_x_band = lim_otf_x[:num_keep]
    lim_otf_x_band = np.mean(lim_otf_x_band.reshape(-1, lim_band_width), axis=1)
    lim_otf_y_band = lim_otf_y[:num_keep]
    lim_otf_y_band_avg = np.mean(lim_otf_y_band.reshape(-1, lim_band_width), axis=1)
    lim_otf_y_band_std = np.std(lim_otf_y_band.reshape(-1, lim_band_width), axis=1)

    # ITF/OTF.
    num_keep = lim_band_width*int(len(lim_com_x) / lim_band_width)

# Plot it.
if plot_opt == 1:
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(12, 5))
    ax1.plot(lams_itf_x, lams_itf_y, label='LAMS', color=tableau20[12], lw=3)
    ax1.plot(lim_itf_x,  lim_itf_y, label='3DLIM', color=tableau20[18], lw=3)
    ax1.plot(rbs_itf_x,  rbs_itf_y, '*', label='RBS', color=tableau20[12], mec='k', ms=12)
    ax2.plot(lams_otf_x, lams_otf_y, label='LAMS', color=tableau20[12], lw=3)
    ax2.plot(lim_otf_x,  lim_otf_y, label='3DLIM', color=tableau20[18], lw=3)
    ax2.plot(rbs_otf_x,  rbs_otf_y, '*', label='RBS', color=tableau20[12], mec='k', ms=12)
    if log_plot:
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.set_ylim([0.001, 3])
        ax2.set_ylim([0.001, 3])
    else:
        ax1.set_ylim([0, 1.1])
        ax2.set_ylim([0, 1.1])
    ax1.set_xlabel('Distance along probe (cm)', fontsize=16)
    ax2.set_xlabel('Distance along probe (cm)', fontsize=16)
    ax1.set_ylabel('Deposition (normalized)', fontsize=16)
    ax1.set_xlim([0, 8])
    ax2.legend(fontsize=16)
    ax1.annotate('ITF', (0.5, 0.9), xycoords='axes fraction', fontsize=24)
    ax2.annotate('OTF', (0.5, 0.9), xycoords='axes fraction', fontsize=24)
    ax1.tick_params(labelsize=12)
    ax2.tick_params(labelsize=12)
    fig.tight_layout()
    fig.show()

    # Plot the ITF/OTF ratio along for LAMS and 3DLIM.
    if itfotf_plot:

        fig, ax = plt.subplots()
        ax.plot(lams_com_x, lams_itfotf, color=tableau20[12])
        ax.plot(lim_com_x,  lim_itfotf,  color=tableau20[18])
        ax.set_ylabel("ITF/OTF")
        ax.set_xlabel("Distance along probe (cm)")
        fig.tight_layout()
        fig.show()

elif plot_opt == 2:

    # Color scheme.
    plt.style.use("tableau-colorblind10")

    # Plot constants.
    fontsize  = 14
    labelsize = 11
    band_alpha = 0.6
    lam_color = "C5"
    rbs_color = "C5"
    lim_color = "C4"
    ms = 12
    lw = 3
    fontname = 'Arial'

    # Grid it.
    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(2, 1, 2)
    ax2 = fig.add_subplot(2, 2, 1)
    ax3 = fig.add_subplot(2, 2, 2)

    # Remove top and right spines.
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)

    # Set log scale.
    if log_plot:
        ax2.set_yscale("log")
        ax3.set_yscale("log")
        ax2.set_ylim([0.005, 3])
        ax3.set_ylim([0.005, 3])

    # Plot the LAMS data.
    if lam_band:
        ax2.fill_between(lams_itf_x, lams_itf_y_errp, lams_itf_y_errm, alpha=band_alpha)
        ax3.fill_between(lams_otf_x, lams_otf_y_errp, lams_otf_y_errm, alpha=band_alpha)
    else:
        ax1.plot(lams_com_x, lams_itfotf, color=lam_color, lw=lw)
        ax2.plot(lams_itf_x, lams_itf_y, color=lam_color, lw=lw)
        ax3.plot(lams_otf_x, lams_otf_y, color=lam_color, lw=lw, label="LAMS")

    # Plot the 3DLIM data.
    if lim_band:
        x = lim_itf_x_band; y = lim_itf_y_band_avg; err = lim_itf_y_band_std
        ax2.fill_between(x, y-err, y+err, alpha=band_alpha)
        x = lim_otf_x_band; y = lim_otf_y_band_avg; err = lim_otf_y_band_std
        ax3.fill_between(x, y-err, y+err, alpha=band_alpha)
    else:
        ax1.plot(lim_com_x, lim_itfotf, color=lim_color, lw=lw)
        ax2.plot(lim_itf_x, lim_itf_y, color=lim_color, lw=lw)
        ax3.plot(lim_otf_x, lim_otf_y, color=lim_color, lw=lw, label="3DLIM")

    # Plot the RBS data.
    ax2.plot(rbs_itf_x, rbs_itf_y, '*', color=rbs_color, ms=ms, mec='k')
    ax3.plot(rbs_otf_x, rbs_otf_y, '*', color=rbs_color, ms=ms, mec='k', label="RBS")

    # Set limits.
    ax1.set_xlim([0, 8])
    ax1.set_ylim([0, 5])
    ax2.set_xlim([0, 8])
    ax3.set_xlim([0, 8])

    # Don't need tick labels on the OTF one.
    ax3.tick_params(labelleft=False)

    # Set explicit tick labels.
    ax1.set_xticks(np.arange(0, 10, 2))
    ax2.set_xticks(np.arange(0, 10, 2))
    ax3.set_xticks(np.arange(0, 10, 2))

    # Labels and font adjustments.
    ax1.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
    ax1.set_ylabel("ITF/OTF\n", fontsize=fontsize)
    ax2.set_ylabel("Deposition (normalized)", fontsize=fontsize)
    ax1.tick_params(which='both', labelsize=labelsize)
    ax2.tick_params(which='both', labelsize=labelsize)
    ax3.tick_params(which='both', labelsize=labelsize)

    # Annotations.
    ax2.text(0.15, 0.15, 'ITF', transform=ax2.transAxes, fontsize=18)
    ax3.text(0.15, 0.15, 'OTF', transform=ax3.transAxes, fontsize=18)
    ax3.legend(fontsize=fontsize)

    fig.tight_layout()
    fig.show()
