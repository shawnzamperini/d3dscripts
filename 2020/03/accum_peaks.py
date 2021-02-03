import matplotlib.pyplot as plt
import numpy as np
import oedge_plots
import tkinter as tk
from tkinter import filedialog
from scipy.optimize import curve_fit


# Some constants.
ring = 30
smin = 15
smax = 45

# First load in the OedgePlots object with all the data.
root = tk.Tk(); root.withdraw()
netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
op = oedge_plots.OedgePlots(netcdf_path)
op.add_dat_file(netcdf_path.split('.nc')[0] + '.dat')
print(netcdf_path)

# Grab just the x, y data of the along ring plot of s vs. imp. density.
x, y = op.along_ring(ring, 'DDLIMS', charge='all', plot_it=False)
y = y / 1e14

# Plot it up, asking the user use to zoom in on the impurity peak. Then pull
# those plot limits to grab just the subset of data the will be fit to a peak.
fig, ax = plt.subplots()
ax.plot(x, y, 'k')
ax.set_ylim([0, None])
ax.set_xlabel('S (m)', fontsize=16)
ax.set_ylabel('Impurity Density (m-3)', fontsize=16)
fig.tight_layout()
fig.show()

ans = input('Zoom in on impurity peak and press enter when done...')
xlims = ax.get_xlim()

peak_idx = np.where(np.logical_and(x > xlims[0], x < xlims[1]))
x_peak = x[peak_idx]; y_peak = y[peak_idx]

def gaussian(x, amp1, cen1, sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

#guess = (x_peak.max() / 2, y_peak.max(), 5)
guess = (y_peak.max() / 2, x_peak.max(), x_peak.max() / 4)
popt, pcov = curve_fit(gaussian, x_peak, y_peak, p0=guess)

peak_val = gaussian(popt[1], *popt)
actual_maxy = y_peak.max()
actual_maxx = x_peak[y_peak == actual_maxy]

y_fit = gaussian(x_peak, *popt)
ax.plot(x_peak, y_fit, 'r')
ax.plot(actual_maxx, actual_maxy, 'k.', ms=10)
ax.plot(popt[1], peak_val, 'r.', ms=10)
fig.tight_layout()
fig.show()

# Estimate average impurity content between X-points.
idx = np.where(np.logical_and(x>=smin, x<=smax))[0]
avg_imp = y[idx].mean()

print()
print('Actual Values:')
print('Peak Value: {:.2e}'.format(actual_maxy))
#print('Peak Value: {:.2e}'.format(actual_maxy*1e14))
print('Peak S:     {:.2f}'.format(actual_maxx[0]))
print('Average:    {:.2e}'.format(avg_imp))
print()
print('Gaussian Values:')
print('Peak Value: {:.2e}'.format(peak_val*1e14))
print('Amplitude:  {:.2e}'.format(popt[0]))
print('S Center:   {:.2f}'.format(popt[1]))
print('Width:      {:.2f}'.format(popt[2]))
