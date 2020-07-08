import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pretty_plots as pp
import sys
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


# Constant to go from dist to R-rsep OMP
slope = 1.1648
intercept = 9.9438

# LAMS calibrations
cal_sl = 4.952E-7

cu05_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/CU05_Map_Analysis2.xlsx'
cd05_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/CD05_Map_Analysis2.xlsx'

dfu = pd.read_excel(cu05_file)
dfd = pd.read_excel(cd05_file)

# Do the pivot thing to have the index as axial values.
dfu = dfu.pivot(columns='z Location [mm]', index='Axial Location [mm]')
dfd = dfd.pivot(columns='z Location [mm]', index='Axial Location [mm]')

# Then the LAMS values.
u_x     = dfu.index.values / 10  # mm to cm
#u_y     = dfu.mean(axis=1).values
u_y     = dfu[[('Total W', 0.75), ('Total W', 1.0), ('Total W', 1.25), ('Total W', 1.5), ('Total W', 1.75)]].mean(axis=1).values
#u_y_err = dfu.std(axis=1).values
u_y_err = dfu[[('Total W', 0.75), ('Total W', 1.0), ('Total W', 1.25), ('Total W', 1.5), ('Total W', 1.75)]].std(axis=1).values
d_x     = dfd.index.values / 10  # mm to cm
#d_y     = dfd.mean(axis=1).values
d_y     = dfd[[('Total W', 0.75), ('Total W', 1.0), ('Total W', 1.25), ('Total W', 1.5), ('Total W', 1.75)]].mean(axis=1).values
d_y_err = dfd[[('Total W', 0.75), ('Total W', 1.0), ('Total W', 1.25), ('Total W', 1.5), ('Total W', 1.75)]].std(axis=1).values
#d_y_err = dfd.std(axis=1).values

# Assign o correct side
itf_x = d_x
itf_y = d_y
itf_y_err = d_y_err
otf_x = u_x
otf_y = u_y
otf_y_err = u_y_err

# Do interpolation functions so we can divide the two sides.
f_itf     = interp1d(itf_x, itf_y)
f_itf_err = interp1d(itf_x, itf_y_err)
f_otf     = interp1d(otf_x, otf_y)
f_otf_err = interp1d(otf_x, otf_y_err)

x_int       = np.linspace(0, 5, 500)
itf_int     = f_itf(x_int)
itf_err_int = f_itf_err(x_int)
otf_int     = f_otf(x_int)
otf_err_int = f_otf_err(x_int)

rminrsepomp = slope * x_int + intercept

# Calculate ITF/OTF.
itfotf = itf_int / otf_int
itfotf_err = itfotf * np.sqrt(np.square(itf_err_int/itf_int)
             + np.square(otf_err_int/otf_int))

# Smooth it.
itfotf_filt = savgol_filter(itfotf, 101, 3)
itfotf_err_filt = savgol_filter(itfotf_err, 101, 3)

# Print out values.
"""
print('Distance along probe:')
for i in x_int:
    print(i)
print('ITF/OTF')
for i in itfotf:
    print(i)
print('ITF/OTF Error')
for i in itfotf_err:
    print(i)
"""

# Plot it.
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(10,5))
ax1.plot(rminrsepomp, itf_int*cal_sl, 'r', label='ITF')
ax1.plot(rminrsepomp, otf_int*cal_sl, 'b', label='OTF')
#ax2.plot(rminrsepomp, itfotf)
ax2.plot(rminrsepomp, itfotf_filt)
ax2.fill_between(rminrsepomp, itfotf-itfotf_err, itfotf+itfotf_err, alpha=0.2)
ax1.set_xlabel('R-Rsep OMP (cm)', fontsize=16)
ax2.set_xlabel('R-Rsep OMP (cm)', fontsize=16)
ax1.legend(fontsize=16)
ax1.set_ylabel('W Areal Density (1e15 cm2)', fontsize=16)
ax2.set_ylabel('ITF/OTF', fontsize=16)
#ax2.set_xlim([10, 14])
ax2.set_ylim([0,2])
fig.tight_layout()
fig.show()
