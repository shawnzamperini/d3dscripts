import matplotlib.pyplot as plt
import numpy as np
import lim_plots
from scipy.signal import savgol_filter


# Some constants.
rsepx1 = 1.0907
rsepx2 = 6.9035
drop_tip = 5
smooth = True

# Load a diffusive and convective case.
diff_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/archive/colprobe-z2-060d.nc'
#conv_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/colprobe-z2-061e.nc'
conv_path = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/archive/colprobe-z2-061h.nc'

# Data for A2 since it's common to compare results against it.
a2_itf_x = np.array([8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3,
                     2.5, 2, 1.5, 1, 0.5, 0])
a2_itf_y = np.array([0.001363767, 0.001363767, 0.0040913, 0, 0.002727533,
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

# Normalize the RBS data.
max_a2_y = max(a2_itf_y.max(), a2_otf_y.max())
a2_itf_y = a2_itf_y / max_a2_y
a2_otf_y = a2_otf_y / max_a2_y

# Load centerline 3DLIM data.
lp_diff = lim_plots.LimPlots(diff_path)
lp_conv = lim_plots.LimPlots(conv_path)
diff_cent = lp_diff.centerline(show_plot=False)
conv_cent = lp_conv.centerline(show_plot=False)

# Get the data, dropping the tip values.
otf_x_diff = diff_cent['otf_x'][:-drop_tip]
otf_y_diff = diff_cent['otf_y'][:-drop_tip]
itf_x_diff = diff_cent['itf_x'][drop_tip:]
itf_y_diff = diff_cent['itf_y'][drop_tip:]
otf_x_conv = conv_cent['otf_x'][:-drop_tip]
otf_y_conv = conv_cent['otf_y'][:-drop_tip]
itf_x_conv = conv_cent['itf_x'][drop_tip:]
itf_y_conv = conv_cent['itf_y'][drop_tip:]

# Scale the X values to R-Rsep OMP.
otf_x_diff = rsepx1 * otf_x_diff * 100 + rsepx2
itf_x_diff = rsepx1 * itf_x_diff * 100 + rsepx2
otf_x_conv = rsepx1 * otf_x_conv * 100 + rsepx2
itf_x_conv = rsepx1 * itf_x_conv * 100 + rsepx2
a2_itf_x = rsepx1 * a2_itf_x + rsepx2
a2_otf_x = rsepx1 * a2_otf_x + rsepx2

# Normalize the 3DLIM data.
max_diff = max(otf_y_diff.max(), itf_y_diff.max())
otf_y_diff = otf_y_diff / max_diff
itf_y_diff = itf_y_diff / max_diff
max_conv = max(otf_y_conv.max(), itf_y_conv.max())
otf_y_conv = otf_y_conv / max_conv
itf_y_conv = itf_y_conv / max_conv

# Then scale all the Y data back to areal density.
otf_y_diff *= max_a2_y * 1e15
itf_y_diff *= max_a2_y * 1e15
otf_y_conv *= max_a2_y * 1e15
itf_y_conv *= max_a2_y * 1e15
a2_otf_y *= max_a2_y * 1e15
a2_itf_y *= max_a2_y * 1e15

# Smooth the data.
if smooth:
    window = 21
    otf_y_diff = savgol_filter(otf_y_diff, window, 3)
    itf_y_diff = savgol_filter(itf_y_diff, window, 3)
    otf_y_conv = savgol_filter(otf_y_conv, window, 3)
    itf_y_conv = savgol_filter(itf_y_conv, window, 3)

# Plotting.
color1 = "tab:red"
color2 = "tab:purple"
color3 = "tab:cyan"
fontsize = 16
ms = 17
lw = 3
fig, (ax1, ax2) = plt.subplots(figsize=(10,5), ncols=2, nrows=1)
ax1.plot(otf_x_diff, otf_y_diff, lw=lw+2, color='k')
ax1.plot(otf_x_diff, otf_y_diff, lw=lw, color=color1)
ax2.plot(itf_x_diff, itf_y_diff, lw=lw+2, color='k')
ax2.plot(itf_x_diff, itf_y_diff, lw=lw, label="Diffusive", color=color1)
ax1.plot(otf_x_conv, otf_y_conv, lw=lw+2, color='k')
ax1.plot(otf_x_conv, otf_y_conv, lw=lw, color=color2)
ax2.plot(itf_x_conv, itf_y_conv, lw=lw+2, color='k')
ax2.plot(itf_x_conv, itf_y_conv, lw=lw, label="Convective", color=color2)
#ax1.plot(a2_otf_x, a2_otf_y, '*', ms=ms, mec='k', mew=1, color=color3)
#ax2.plot(a2_itf_x, a2_itf_y, '*', ms=ms, mec='k', mew=1, label="RBS", color=color3)
ax1.errorbar(a2_otf_x, a2_otf_y, yerr=a2_otf_y*0.25, marker='*', ms=ms, mec='k', mew=1.5, lw=0, elinewidth=2, capsize=3, ecolor='k', color=color3)
ax2.errorbar(a2_itf_x, a2_itf_y, yerr=a2_itf_y*0.25, marker='*', ms=ms, mec='k', mew=1.5, lw=0, elinewidth=2, capsize=3, ecolor='k', label="RBS", color=color3)
ax1.set_xlim([6.5, 13])
ax2.set_xlim([6.5, 13])
ax1.set_ylim([1e12, 1e15])
ax2.set_ylim([1e12, 1e15])
ax1.set_yscale("log")
ax2.set_yscale("log")
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.tick_params(labelleft=False)
ax2.legend(loc="lower left", fontsize=14)
ax1.set_xlabel("R-Rsep OMP (cm)", fontsize=14)
ax2.set_xlabel("R-Rsep OMP (cm)", fontsize=14)
ax1.set_ylabel(r"W Areal Density $\mathrm{(W/cm^2)}$", fontsize=fontsize)
fig.tight_layout()
fig.show()

# Save the figure.
fig.savefig('conv_vs_diff.pdf', format='pdf')
