import lim_plots
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import statsmodels.api as sm


# Some constants.
cal = 0.5E-06  # LAMS calibration to go from counts to areal density.
tip_itf_ignore = 0  # Ignore the first this many points on 3DLIM.
tip_otf_ignore = 0
tip_itf_lams_ignore = 0  # Ignore the first this many point on LAMS.
tip_otf_lams_ignore = 0
#log_plot = False
log_plot = True
show_unused = True
itfotf_plot = True
plot_opt = 3
lim_band = False
lim_band_width = 5
lam_band = False
smooth         = True
smooth_lams    = True
smooth_type    = "savgol"
#smooth_type    = "lowess"
frac           = 1/6      # For lowess smoothing.
window         = 21       # For savgol filtering.
window_lams    = 41
to_omp         = True

# Path to the 3DLIM data.
ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a15-001g.nc'

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
a15_rbs = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
         'Collector Probe Excel Sheets/A15.xlsx'

# Path to A8 LAMS data.
au15_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AU15_Map_Analysis.xlsx'
ad15_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AD15_Map_Analysis.xlsx'

# Load 3DLIM data.
lp = lim_plots.LimPlots(ncpath)

# Load into Dataframes.
rbs = pd.read_excel(a15_rbs, nrows=20)
au15_lams = pd.read_excel(au15_file)
ad15_lams = pd.read_excel(ad15_file)

# Pivot the data so the index is the poloidal data and the columns the radial.
au15_lams = au15_lams.pivot(columns='z Location [mm]', index='Axial Location [mm]')
ad15_lams = ad15_lams.pivot(columns='z Location [mm]', index='Axial Location [mm]')

# Pull out the data we will plot. RBS first.
rbs_itf_x     = rbs['Distance from Tip U (cm)'].values
rbs_itf_y     = rbs['W Areal Density U (1e15 W/cm2)'].values
rbs_itf_y_err = rbs['W Areal Density Error U (1e15 W/cm2)'].values
rbs_otf_x     = rbs['Distance from Tip D (cm)'].values
rbs_otf_y     = rbs['W Areal Density D (1e15 W/cm2)'].values
rbs_otf_y_err = rbs['W Areal Density Error D (1e15 W/cm2)'].values

# Then the LAMS values. Forward: U=ITF, D=OTF.
#lams_itf_x     = au15_lams.index.values / 10  # mm to cm
#lams_itf_y     = au15_lams.mean(axis=1).values
#lams_itf_y_err = au15_lams.std(axis=1).values
#lams_otf_x     = ad15_lams.index.values / 10  # mm to cm
#lams_otf_y     = ad15_lams.mean(axis=1).values
#lams_otf_y_err = ad15_lams.std(axis=1).values

# Get the centerline data.
lams_itf_x = au15_lams["Total W"][2.5].index.values / 10  # mm to cm
lams_itf_y = au15_lams["Total W"][2.5].values
lams_itf_y_err = au15_lams["Total W error"][2.5].values
lams_otf_x = ad15_lams["Total W"][2.5].index.values / 10  # mm to cm
lams_otf_y = ad15_lams["Total W"][2.5].values
lams_otf_y_err = ad15_lams["Total W error"][2.5].values

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

# Drop some RBS points we don't want.
rbs_otf_drop = [0, 3, 4, 8, 15, 16]
rbs_otf_x = np.delete(rbs_otf_x, rbs_otf_drop)
rbs_otf_y = np.delete(rbs_otf_y, rbs_otf_drop)
rbs_otf_y_err = np.delete(rbs_otf_y_err, rbs_otf_drop)
rbs_itf_drop = [1,2,3,5,6,7,8,10,12,16,17,18]
rbs_itf_x = np.delete(rbs_itf_x, rbs_itf_drop)
rbs_itf_y = np.delete(rbs_itf_y, rbs_itf_drop)
rbs_itf_y_err = np.delete(rbs_itf_y_err, rbs_itf_drop)

# On the OTF side linearly interpolated across the scraped region.
scraped = np.where(np.logical_and(lams_otf_x>=1.3, lams_otf_x<=2.7))[0]
x1 = lams_otf_x[scraped[0]]
x2 = lams_otf_x[scraped[-1]]
y1 = lams_otf_y[scraped[0]]
y2 = lams_otf_y[scraped[-1]]
m = (y2 - y1) / (x2 - x1)
x_int = np.linspace(x1, x2, len(scraped))
y_int = m * (x_int - x1) + y1
lams_otf_y[scraped] = y_int

# For A15, need to restrict the data to everything after 1.5 cm due to scraping :(
start = 1.5
rbs_itf_keep = np.where(rbs_itf_x>=start)[0]
rbs_otf_keep = np.where(rbs_otf_x>=start)[0]
lams_itf_keep = np.where(lams_itf_x>=start)[0]
lams_otf_keep = np.where(lams_otf_x>=start)[0]
lim_itf_keep = np.where(lim_itf_x>=start)[0]
lim_otf_keep = np.where(lim_otf_x>=start)[0]
rbs_itf_x = rbs_itf_x[rbs_itf_keep]
rbs_itf_y = rbs_itf_y[rbs_itf_keep]
rbs_itf_y_err = rbs_itf_y_err[rbs_itf_keep]
rbs_otf_x = rbs_otf_x[rbs_otf_keep]
rbs_otf_y = rbs_otf_y[rbs_otf_keep]
rbs_otf_y_err = rbs_otf_y_err[rbs_otf_keep]

# Save the not considered LAMS and 3DLIM data.
lams_itf_x_unused = np.delete(lams_itf_x, lams_itf_keep)
lams_itf_y_unused = np.delete(lams_itf_y, lams_itf_keep)
lams_itf_y_err_unused = np.delete(lams_itf_y_err, lams_itf_keep)
lams_otf_x_unused = np.delete(lams_otf_x, lams_otf_keep)
lams_otf_y_unused = np.delete(lams_otf_y, lams_otf_keep)
lams_otf_y_err_unused = np.delete(lams_otf_y_err, lams_otf_keep)
lim_itf_x_unused = np.delete(lim_itf_x, lim_itf_keep)
lim_itf_y_unused = np.delete(lim_itf_y, lim_itf_keep)
#lim_itf_y_err_unused = np.delete(lim_itf_y_err, lim_itf_keep)
lim_otf_x_unused = np.delete(lim_otf_x, lim_otf_keep)
lim_otf_y_unused = np.delete(lim_otf_y, lim_otf_keep)
#lim_otf_y_err_unused = np.delete(lim_otf_y_err, lim_otf_keep)

# The rest of restricting the data.
lams_itf_x = lams_itf_x[lams_itf_keep]
lams_itf_y = lams_itf_y[lams_itf_keep]
lams_itf_y_err = lams_itf_y_err[lams_itf_keep]
lams_otf_x = lams_otf_x[lams_otf_keep]
lams_otf_y = lams_otf_y[lams_otf_keep]
lams_otf_y_err = lams_otf_y_err[lams_otf_keep]
lim_itf_x = lim_itf_x[lim_itf_keep]
lim_itf_y = lim_itf_y[lim_itf_keep]
lim_otf_x = lim_otf_x[lim_otf_keep]
lim_otf_y = lim_otf_y[lim_otf_keep]

# Ignore number of data points.
lim_itf_x = lim_itf_x[tip_itf_ignore:]
lim_itf_y = lim_itf_y[tip_itf_ignore:]
lim_otf_x = lim_otf_x[tip_otf_ignore:]
lim_otf_y = lim_otf_y[tip_otf_ignore:]
lams_itf_x = lams_itf_x[tip_itf_lams_ignore:]
lams_itf_y = lams_itf_y[tip_itf_lams_ignore:]
lams_otf_x = lams_otf_x[tip_otf_lams_ignore:]
lams_otf_y = lams_otf_y[tip_otf_lams_ignore:]

# Remove spurious spikes in the data. Best to index from the end since the
# tip_ignore constants can change these values when counting from the front,
# but not from the back.
#lams_itf_y[23] = (lams_itf_y[22] + lams_itf_y[24]) / 2
#lams_itf_y[-377] = (lams_itf_y[-378] + lams_itf_y[-376]) / 2
#lams_otf_y[165] = (lams_otf_y[164] + lams_otf_y[166]) / 2
lams_itf_y[-174] = (lams_itf_y[-175] + lams_itf_y[-173]) / 2

# Linearly interpolated across a spiked region.
#x1 = lams_itf_x[130]
#x2 = lams_itf_x[144]
#y1 = lams_itf_y[130]
#y2 = lams_itf_y[144]
#x1 = lams_otf_x[-270]
#x2 = lams_otf_x[-256]
#y1 = lams_otf_y[-270]
#y2 = lams_otf_y[-256]
#m = (y2 - y1) / (x2 - x1)
#x_int = np.linspace(x1, x2, 13)
#y_int = m * (x_int - x1) + y1
#lams_itf_y[131:144] = y_int
#lams_otf_y[-269:-256] = y_int

# Normalize all the data.
max_rbs = max(rbs_itf_y.max(), rbs_otf_y.max())
rbs_itf_y = rbs_itf_y / max_rbs
rbs_otf_y = rbs_otf_y / max_rbs
max_lams = max(lams_itf_y.max(), lams_otf_y.max())
lams_itf_y = lams_itf_y / max_lams
lams_otf_y = lams_otf_y / max_lams
max_lim = max(lim_itf_y.max(), lim_otf_y.max())
lim_itf_y = lim_itf_y / max_lim
lim_otf_y = lim_otf_y / max_lim

# Including the unused data.
lams_itf_y_unused = lams_itf_y_unused / max_lams
lams_otf_y_unused = lams_otf_y_unused / max_lams
lim_itf_y_unused = lim_itf_y_unused / max_lim
lim_otf_y_unused = lim_otf_y_unused / max_lim

# Can append the lim data back onto our array now.
#lim_otf_x = np.append(lim_otf_x_unused, lim_otf_x)
#lim_otf_y = np.append(lim_otf_y_unused, lim_otf_y)
#lim_itf_x = np.append(lim_itf_x_unused, lim_itf_x)
#lim_itf_y = np.append(lim_itf_y_unused, lim_itf_y)
lim_itf_x_all = np.append(lim_itf_x_unused, lim_itf_x)
lim_otf_x_all = np.append(lim_otf_x_unused, lim_otf_x)
lim_itf_y_all = np.append(lim_itf_y_unused, lim_itf_y)
lim_otf_y_all = np.append(lim_otf_y_unused, lim_otf_y)

# Smooth out the data.
if smooth:
    if smooth_type == "savgol":
        lim_otf_y_all = savgol_filter(lim_otf_y_all, window, 3)
        lim_itf_y_all = savgol_filter(lim_itf_y_all, window, 3)
    elif smooth_type == "lowess":
        lowess = sm.nonparametric.lowess
        lim_otf_y_all = lowess(lim_otf_y_all, lim_otf_x_all, return_sorted=False, frac=frac)
        lim_itf_y_all = lowess(lim_itf_y_all, lim_itf_x_all, return_sorted=False, frac=frac)
    else:
        print("Error: Incorrect smooth_type.")

    lams_otf_y = savgol_filter(lams_otf_y, window_lams, 3)
    lams_itf_y = savgol_filter(lams_itf_y, window_lams, 3)

# ITF/OTF plot: Need to interpolate onto common x values to calculate ratios.
#lim_itf_x_all = np.append(lim_itf_x_unused, lim_itf_x)
#lim_otf_x_all = np.append(lim_otf_x_unused, lim_otf_x)
#lim_itf_y_all = np.append(lim_itf_y_unused, lim_itf_y)
#lim_otf_y_all = np.append(lim_otf_y_unused, lim_otf_y)
num_points = len(lams_itf_x)
#num_points = 300
lams_com_x = np.linspace(max(lams_itf_x.min(), lams_otf_x.min()), min(lams_itf_x.max(), lams_otf_x.max()), num_points)
lim_com_x  = np.linspace(max(lim_itf_x_all.min(),  lim_otf_x_all.min()),  min(lim_itf_x_all.max(),  lim_otf_x_all.max()),  num_points+len(lim_itf_x_unused))
f_lams_itf = interp1d(lams_itf_x, lams_itf_y)
f_lams_otf = interp1d(lams_otf_x, lams_otf_y)
f_lim_itf  = interp1d(lim_itf_x_all,  lim_itf_y_all)
f_lim_otf  = interp1d(lim_otf_x_all,  lim_otf_y_all)
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

# Convert to R-Rsep OMP coordinates.
if to_omp:
    m = 1.0907
    y1 = 6.9035
    lams_com_x = m * lams_com_x + y1
    lim_com_x  = m * lim_com_x  + y1
    lams_itf_x = m * lams_itf_x + y1
    lams_otf_x = m * lams_otf_x + y1
    lim_itf_x  = m * lim_itf_x  + y1
    lim_otf_x  = m * lim_otf_x  + y1
    rbs_itf_x  = m * rbs_itf_x  + y1
    rbs_otf_x  = m * rbs_otf_x  + y1
    lim_itf_x_all = m * lim_itf_x_all + y1
    lim_otf_x_all = m * lim_otf_x_all + y1
    xlabel = "R-Rsep OMP (cm)"
    xbounds = [6.5, 15.5]
    xlabels = np.arange(7, 17, 2)

else:
    xlabel = "Distance along probe (cm)"
    xbounds = [0, 8]
    xlabels = np.arange(0, 10, 2)

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
        ax2.fill_between(lams_itf_x, lams_itf_y-lams_itf_y_err, lams_itf_y+lams_itf_y_err, alpha=band_alpha)
        ax3.fill_between(lams_otf_x, lams_otf_y-lams_otf_y_err, lams_otf_y+lams_otf_y_err, alpha=band_alpha)
    else:
        ax1.plot(lams_com_x, lams_itfotf, color=lam_color, lw=lw)
        ax2.plot(lams_itf_x, lams_itf_y, color=lam_color, lw=lw)
        ax3.plot(lams_otf_x, lams_otf_y, color=lam_color, lw=lw, label="LAMS")
        ax2.plot(lams_itf_x_unused, lams_itf_y_unused, color=lam_color, lw=lw, linestyle='--')
        ax3.plot(lams_otf_x_unused, lams_otf_y_unused, color=lam_color, lw=lw, linestyle='--')

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
        ax2.plot(lim_itf_x_unused, lim_itf_y_unused, color=lim_color, lw=lw)
        ax3.plot(lim_otf_x_unused, lim_otf_y_unused, color=lim_color, lw=lw)

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

elif plot_opt == 3:

    # Color scheme.
    plt.style.use("tableau-colorblind10")

    # Plot constants.
    fontsize  = 14
    labelsize = 11
    band_alpha = 0.6
    lam_color = "C5"
    rbs_color = "C5"
    lim_color = "C4"
    lam_color = "lightskyblue"
    rbs_color = "lightskyblue"
    #lim_color = "tab:purple"
    lim_color = "lightskyblue"
    ms = 12
    lw = 4
    fontname = 'Arial'

    # Grid it.
    fig = plt.figure(figsize=(10, 4))
    ax2 = fig.add_subplot(1, 2, 1)
    ax3 = fig.add_subplot(1, 2, 2)

    # Remove top and right spines.
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)

    # Set log scale.
    if log_plot:
        ax2.set_yscale("log")
        ax3.set_yscale("log")
        ax2.set_ylim([0.005*max_rbs, 3*max_rbs])
        ax3.set_ylim([0.005*max_rbs, 3*max_rbs])

    # Plot the LAMS data.
    if lam_band:
        #ax2.fill_between(lams_itf_x, lams_itf_y_errp, lams_itf_y_errm, alpha=band_alpha, color=lam_color)
        #ax3.fill_between(lams_otf_x, lams_otf_y_errp, lams_otf_y_errm, alpha=band_alpha, color=lam_color)
        ax2.fill_between(lams_itf_x, (lams_itf_y+lams_itf_y*0.2)*max_rbs, (lams_itf_y-lams_itf_y*0.2)*max_rbs, alpha=band_alpha, color=lam_color)
        ax3.fill_between(lams_otf_x, (lams_otf_y+lams_otf_y*0.2)*max_rbs, (lams_otf_y-lams_otf_y*0.2)*max_rbs, alpha=band_alpha, color=lam_color, label="LAMS")
    else:
        ax2.plot(lams_itf_x, lams_itf_y*max_rbs, color=lam_color, lw=lw)
        ax3.plot(lams_otf_x, lams_otf_y*max_rbs, color=lam_color, lw=lw, label="LAMS")

    # Plot the 3DLIM data.
    if lim_band:
        x = lim_itf_x_band; y = lim_itf_y_band_avg; err = lim_itf_y_band_std
        ax2.fill_between(x, y-err, y+err, alpha=band_alpha)
        x = lim_otf_x_band; y = lim_otf_y_band_avg; err = lim_otf_y_band_std
        ax3.fill_between(x, y-err, y+err, alpha=band_alpha)
    else:
        ax2.plot(lim_itf_x_all, lim_itf_y_all*max_rbs, color='k', lw=lw+2)
        ax2.plot(lim_itf_x_all, lim_itf_y_all*max_rbs, color=lim_color, lw=lw)
        ax3.plot(lim_otf_x_all, lim_otf_y_all*max_rbs, color='k', lw=lw+2)
        ax3.plot(lim_otf_x_all, lim_otf_y_all*max_rbs, color=lim_color, lw=lw, label="3DLIM")

    # Plot the RBS data.
    ax2.plot(rbs_itf_x, rbs_itf_y*max_rbs, '*', color=rbs_color, ms=ms, mec='k')
    ax3.plot(rbs_otf_x, rbs_otf_y*max_rbs, '*', color=rbs_color, ms=ms, mec='k', label="RBS")

    # Set limits.
    #ax2.set_xlim([0, 8])
    #ax3.set_xlim([0, 8])
    ax2.set_xlim(xbounds)
    ax3.set_xlim(xbounds)

    # Don't need tick labels on the OTF one.
    ax3.tick_params(labelleft=False)

    # Set explicit tick labels.
    #ax2.set_xticks(np.arange(0, 10, 2))
    #ax3.set_xticks(np.arange(0, 10, 2))
    ax2.set_xticks(xlabels)
    ax3.set_xticks(xlabels)

    # Labels and font adjustments.
    ax2.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
    ax3.set_xlabel("Distance along probe (cm)", fontsize=fontsize)
    ax2.set_ylabel(r"Deposition (1e15 W/$\mathrm{cm^2}$)", fontsize=fontsize)
    ax2.tick_params(which='both', labelsize=labelsize)
    ax3.tick_params(which='both', labelsize=labelsize)

    # Annotations.
    ax2.text(0.15, 0.15, 'ITF', transform=ax2.transAxes, fontsize=18)
    ax3.text(0.15, 0.15, 'OTF', transform=ax3.transAxes, fontsize=18)
    #ax3.legend(fontsize=fontsize)

    fig.tight_layout()
    fig.show()

"""
# Plot it.
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
if show_unused:
    ax1.plot(lams_itf_x_unused, lams_itf_y_unused, color=tableau20[12], lw=3, linestyle=":")
    ax2.plot(lams_otf_x_unused, lams_otf_y_unused, color=tableau20[12], lw=3, linestyle=":")
    ax1.plot(lim_itf_x_unused, lim_itf_y_unused, color=tableau20[18], lw=3)
    ax2.plot(lim_otf_x_unused, lim_otf_y_unused, color=tableau20[18], lw=3)

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

    # Need to interpolate onto common x values to calculate ratios.
    lams_com_x = np.linspace(max(lams_itf_x.min(), lams_otf_x.min()), min(lams_itf_x.max(), lams_otf_x.max()), 300)
    lim_com_x  = np.linspace(max(lim_itf_x.min(),  lim_otf_x.min()),  min(lim_itf_x.max(),  lim_otf_x.max()),  300)
    f_lams_itf = interp1d(lams_itf_x, lams_itf_y)
    f_lams_otf = interp1d(lams_otf_x, lams_otf_y)
    f_lim_itf  = interp1d(lim_itf_x,  lim_itf_y)
    f_lim_otf  = interp1d(lim_otf_x,  lim_otf_y)
    lams_itfotf = f_lams_itf(lams_com_x) / f_lams_otf(lams_com_x)
    lim_itfotf  = f_lim_itf(lim_com_x)   / f_lim_otf(lim_com_x)

    fig, ax = plt.subplots()
    ax.plot(lams_com_x, lams_itfotf, color=tableau20[12])
    ax.plot(lim_com_x,  lim_itfotf,  color=tableau20[18])
    ax.set_ylabel("ITF/OTF")
    ax.set_xlabel("Distance along probe (cm)")
    fig.tight_layout()
    fig.show()
"""
