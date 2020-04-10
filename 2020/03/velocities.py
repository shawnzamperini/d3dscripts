import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog


# Some constants.
ring = 40
charge = 'all'
#charge=10
mz_amu = 183.84
md_amu = 2.01
mz = mz_amu * 1.66e-27  # Mass of W in kg
md = md_amu * 1.66e-27   # Mass of D in kg
col_log = 15

# Option to grab additional reruns and average them together, if available.
mult_runs = False

# First load in the OedgePlots object(s) with all the data.
if mult_runs:

    # Load each OedgePlots object for each run (with their .dat files) and put into list.
    ops = []
    runs = input('Enter series of runs (i.e. "073h" would get all the 073h runs and average the results): ')
    base = '/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z1-167196-'
    for i in range(0, 10):
        netcdf_path = base + runs + str(i) + '.nc'
        try:
            op = oedge_plots.OedgePlots(netcdf_path)
            op.add_dat_file(netcdf_path.split('.nc')[0] + '.dat')
            ops.append(op)
            print(netcdf_path)
        except:
            pass
else:

    # Just select the single DIVIMP run we want.
    root = tk.Tk(); root.withdraw()
    netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
    op = oedge_plots.OedgePlots(netcdf_path)
    op.add_dat_file(netcdf_path.split('.nc')[0] + '.dat')
    print(netcdf_path)
    ops = [op]


# Special option if we want just the average of all the charge states.
if charge == 'all':
    charges = np.arange(1, 74+1)
else:
    charges = [charge]

# Along ring quantities that don't depend on charge so get them now.
s, ti  = ops[0].along_ring(ring, 'KTIBS', plot_it=False)
s, ne  = ops[0].along_ring(ring, 'KNBS',  plot_it=False)
s, vi  = ops[0].along_ring(ring, 'Velocity', plot_it=False)

#nz_all = np.zeros((len(ops), len(s)))
#for j in range(0, len(ops)):
#    s, nz = ops[j].along_ring(ring, 'DDLIMS', plot_it=False, charge='all')
#    nz_all[j] = nz
#nz_avg = np.mean(nz_all, axis=0)

# Array to hold sum of all vti values for each charge state. To be averaged after.
vti_all    = np.zeros((len(charges), len(s)))
velavg_all = np.zeros((len(charges)*len(ops), len(s)))
gamma_all  = np.zeros((len(charges)*len(ops), len(s)))
nz_all     = np.zeros((len(charges)*len(ops), len(s)))

for i in range(0, len(charges)):
    charge = charges[i]

    # Get forces.
    #s, ff = ops[0].along_ring(ring, 'ff', plot_it=False, charge=charge)
    s, fig = ops[0].along_ring(ring, 'fig',   plot_it=False, charge=charge)

    # Slowing down time.
    tau_s = 1.47E13 * mz_amu * ti * np.sqrt(ti / md_amu) / \
            ((1 + md_amu / mz_amu) * ne * np.power(charge, 2) * col_log)

    # vTi is FiG * tau_s / mz.
    vti = fig * tau_s / mz

    for j in range(0, len(ops)):

        # The average velocity for this W charge state on this ring.
        s, velavg = ops[j].along_ring(ring, 'VELAVG', plot_it=False, charge=charge)
        s, nz     = ops[j].along_ring(ring, 'DDLIMS', plot_it=False, charge=charge)
        velavg_all[i * len(ops) + j] = velavg
        nz_all[i * len(ops) + j]     = nz
        gamma_all[i * len(ops) + j]  = velavg * nz

    # Add to array with all values.
    vti_all[i]    = vti
    #velavg_all[i] = velavg

# Get average, ignoring zeros. Don't need to do for vti.
velavg_all[velavg_all==0] = np.nan
nz_all[nz_all==0] = np.nan
gamma_all[gamma_all==0] = np.nan
velavg_final = np.nanmean(velavg_all, axis=0)
nz_final = np.nanmean(nz_all, axis=0)
gamma_final = np.nanmean(gamma_all, axis=0)
vti_final = vti_all.mean(axis=0)

# Calculate the flux of impurities.
#gamma_z = velavg_final * nz_avg

# Plotting commands.
fig, axs = plt.subplots(1,2, figsize=(10,5))

# Velocity plot.
axs[0].plot(s, velavg_final, 'k', label='v')
axs[0].plot(s, vi+vti_final, 'r', label='vi+vTi')
axs[0].plot(s, vi, 'r-.', label='vi')
axs[0].plot(s, vti, 'r--', label='vti')
axs[0].axhline(0, color='k', linestyle='-', alpha=0.3)
axs[0].axvline(s.max()/2, color='k', linestyle='-', alpha=0.3)
axs[0].legend()
axs[0].set_xlabel('S (m)', fontsize=16)
axs[0].set_ylabel('Velocity (m/s)', fontsize=16)

# Gamma plot.
axs[1].plot(s, gamma_final, 'k', label='gamma')
#axs[1].plot(s, vi+vti_final, 'r', label='vi+vTi')
#axs[1].plot(s, vi, 'r-.', label='vi')
#axs[1].plot(s, vti, 'r--', label='vti')
axs[1].axhline(0, color='k', linestyle='-', alpha=0.3)
axs[1].axvline(s.max()/2, color='k', linestyle='-', alpha=0.3)
axs[1].legend()
axs[1].set_xlabel('S (m)', fontsize=16)
axs[1].set_ylabel('Gamma W (m-2 s-1)', fontsize=16)
axs[1].set_ylim([-2e17, 2e17])

fig.tight_layout()
fig.show()
