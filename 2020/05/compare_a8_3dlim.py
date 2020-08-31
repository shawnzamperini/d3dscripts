import lim_plots
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d


# Some constants.
cal = 0.5E-06  # LAMS calibration to go from counts to areal density.
tip_itf_ignore = 3  # Ignore the first this many points on 3DLIM.
tip_otf_ignore = 3
tip_itf_lams_ignore = 0  # Ignore the first this many point on LAMS.
tip_otf_lams_ignore = 0
#log_plot = False
log_plot = True
itfotf_plot = True

# Path to the 3DLIM data.
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-045h.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-045e.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-001v.nc'
ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-001w.nc'
#ncpath = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-a8-002a.nc'

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

# Then the LAMS values.
lams_itf_x     = ad8_lams.index.values / 10  # mm to cm
lams_itf_y     = ad8_lams.mean(axis=1).values
lams_itf_y_err = ad8_lams.std(axis=1).values
lams_otf_x     = au8_lams.index.values / 10  # mm to cm
lams_otf_y     = au8_lams.mean(axis=1).values
lams_otf_y_err = au8_lams.std(axis=1).values

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

# This is just messing around to remove a spurious data point that messes with
# the norm in a8-001v. Shouldn't be an issue with enough statistics.
lim_itf_y[14] = (lim_itf_y[13] + lim_itf_y[15]) / 2

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
