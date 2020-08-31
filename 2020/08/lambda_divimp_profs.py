import matplotlib.pyplot as plt
import numpy as np
import oedge_plots
import tkinter as tk
from tkinter import filedialog
from scipy.optimize import curve_fit


# Some constants.
ring = 65

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

ans = input('Zoom in on exponential region and press enter when done...')
xlims = ax.get_xlim()

# Select only the exponential region.
exp_idx = np.where(np.logical_and(x > xlims[0], x < xlims[1]))
x_exp = x[exp_idx]; y_exp = y[exp_idx]

# Set the x values to start at zero so the fit can work fine.
shift = x_exp.min()
x_exp = x_exp - shift

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

popt, pcov = curve_fit(exp_fit, x_exp, y_exp)

x_fit = np.linspace(x_exp.min(), x_exp.max(), 100)
y_fit = exp_fit(x_fit, *popt)

ax.plot(x_fit + shift, y_fit)
fig.show()

omp_dist = op.nc.variables["MIDIST"][:].data[1][ring]

print("Ring: {}".format(ring))
print("Lambda = {:.2f} cm".format(1/popt[1] * 100))
print("OMP Dist = {:.3f} cm".format(omp_dist * 100))
print("Conn. Length = {:.2f}".format((x_exp+shift).max()))
