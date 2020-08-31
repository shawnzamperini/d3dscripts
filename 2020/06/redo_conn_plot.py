import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# These are the "Tableau 20" colors as RGB. The last one is just black. I added it.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (0, 0, 0)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

xl_path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/07/lp_with_a2.xlsx'
df = pd.read_excel(xl_path, skiprows=1, sheet_name='Sheet1')

all_lp_x = df['R-Rsep OMP Shifted All (cm)'].values[:250]
all_lp_y = df['ne (1e18 m-3)'].values[:250]
bin_lp_x = df['R-Rsep OMP Shifted (cm)'].values[:250]
bin_lp_x_err = df['Error (cm)'].values[:250]
bin_lp_y = df['ne (1e18 m-3).1'].values[:250]
bin_lp_y_err = df['Error'].values[:250]
bin_lp_x = bin_lp_x[~np.isnan(bin_lp_x)][2:]
bin_lp_y = bin_lp_y[~np.isnan(bin_lp_y)][2:]
bin_lp_x_err = bin_lp_x_err[~np.isnan(bin_lp_x_err)][2:]
bin_lp_y_err = bin_lp_y_err[~np.isnan(bin_lp_y_err)][2:]

main_exp_a = 15.796
main_exp_b = -0.239
wind_exp_a = 31441
wind_exp_b = -0.979

main_exp_x = np.linspace(4, 10.3, 100)
wind_exp_x = np.linspace(10.3, 12.5, 100)
main_exp_y = main_exp_a * np.exp(main_exp_b * main_exp_x)
wind_exp_y = wind_exp_a * np.exp(wind_exp_b * wind_exp_x)

xl_path2 = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/Connection Lengths/167196/167196.xlsx'
df2 = pd.read_excel(xl_path2, sheet_name='MAFOT ITF', skiprows=2)
conn_x = df2['R-Rsep OMP (cm)'].values
conn_y = df2['Connection Length (km)'].values * 1000
df3 = pd.read_excel(xl_path2, sheet_name='MAFOT OTF', skiprows=2)
conn_x_otf = df3['R-Rsep OMP (cm)'].values
conn_y_otf = df3['Connection Length (km)'].values * 1000

xl_path3 = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx'
df4 = pd.read_excel(xl_path3)
shift = 1.5
a2_x = df4['R-Rsep omp D (cm)'].values[10:] + shift
a2_y = df4['W Areal Density D (1e15 W/cm2)'].values[10:]
a2_x_err = df4['R-Rsep omp Error D (cm)'].values[10:]
a2_y_err = df4['W Areal Density Error D (1e15 W/cm2)'].values[10:]
a2_x_otf = df4['R-Rsep omp U (cm)'].values[10:] + shift
a2_y_otf = df4['W Areal Density U (1e15 W/cm2)'].values[10:]
a2_x_err_otf = df4['R-Rsep omp Error U (cm)'].values[10:]
a2_y_err_otf = df4['W Areal Density Error U (1e15 W/cm2)'].values[10:]

a2_main_exp_a = 2.708
a2_main_exp_b = -0.201
a2_wind_exp_a = 7.184E7
a2_wind_exp_b = -1.719

a2_main_exp_x = np.linspace(7, 9.5, 100) + shift
a2_main_exp_y = a2_main_exp_a * np.exp(a2_main_exp_b * a2_main_exp_x)
a2_wind_exp_x = np.linspace(10, 12, 100) + shift
a2_wind_exp_y = a2_wind_exp_a * np.exp(a2_wind_exp_b * a2_wind_exp_x)

# Two side by side plots, one with the LP dat and the other with A2.
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(11, 5))
ax1b = ax1.twinx()
ax1.semilogy(main_exp_x, main_exp_y, '--', color='k')
ax1.semilogy(wind_exp_x, wind_exp_y, '--', color='k')
ax1.semilogy(all_lp_x, all_lp_y, '.', alpha=0.1, ms=5, color="k")
ax1.errorbar(bin_lp_x, bin_lp_y, xerr=bin_lp_x_err, yerr=bin_lp_y_err, fmt='k.', ms=10, capsize=3, capthick=1, elinewidth=1, ecolor='k')
plt.yscale("log")
ax1b.semilogy(conn_x, conn_y, color=tableau20[18], lw=3, label="ITF")
ax1b.semilogy(conn_x_otf, conn_y_otf, '--', color=tableau20[18], lw=3, label="OTF")
ax2b = ax2.twinx()
ax2b.semilogy(conn_x_otf, conn_y_otf, '--', color=tableau20[18], lw=3)
ax2b.semilogy(conn_x, conn_y, color=tableau20[18], lw=3)
ax2.semilogy(a2_main_exp_x, a2_main_exp_y, '--', color='k')
ax2.semilogy(a2_wind_exp_x, a2_wind_exp_y, '--', color='k')
ax2.errorbar(a2_x, a2_y, xerr=a2_x_err, yerr=a2_y_err, fmt='k.', ms=12, capsize=3, capthick=1, elinewidth=1, ecolor='k', label="ITF")
ax2.errorbar(a2_x_otf[:-2], a2_y_otf[:-2], xerr=a2_x_err_otf[:-2], yerr=a2_y_err_otf[:-2], fmt='.', mec="k", color="w", ms=12, capsize=3, capthick=1, elinewidth=1, ecolor='k', label="OTF")
plt.yscale("log")
ax1.set_xlim([6, 14])
ax1.set_zorder(10)
ax2.set_zorder(10)
ax1.patch.set_visible(False)
ax2.patch.set_visible(False)
ax1.set_ylabel(r"$\mathrm{n_e}$ ($\mathrm{10^{18}\ m^{-3}}$)", fontsize=16)
ax2.set_ylabel(r"W Areal Density ($\mathrm{10^{15}\ W/cm^{2}}$)", fontsize=16)
ax2b.set_ylabel("Connection Length (m)", fontsize=16, color=tableau20[18])
ax1.tick_params(labelsize=12)
ax2.tick_params(labelsize=12)
ax1b.tick_params(axis='both', which='both', labelsize=12, color=tableau20[18], labelcolor=tableau20[18])
ax2b.tick_params(axis='both', which='both', labelsize=12, color=tableau20[18], labelcolor=tableau20[18])
ax1b.set_ylim([1e-1, 1e2])
ax2b.set_ylim([1e-1, 1e2])
ax1.set_ylim([1e-2, 1e1])
ax2.set_ylim([1e-3, 1e0])
ax2.legend(loc="lower left", fontsize=14)
ax1b.legend(loc="lower left", fontsize=14)
ax1.set_xlabel("R-Rsep OMP (cm)", fontsize=16)
ax2.set_xlabel("R-Rsep OMP (cm)", fontsize=16)
ax1.set_title("Langmuir Probe", fontsize=16)
ax2.set_title("Collector Probe", fontsize=16)
fig.tight_layout()
fig.show()
