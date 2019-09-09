import pretty_plots as pp
import pandas as pd
import numpy as np


xl_path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/07/lp_with_a2.xlsx'
df = pd.read_excel(xl_path, skiprows=1, sheet_name='Sheet1')

all_lp_x = df['R-Rsep OMP (cm)'].values[:250]
all_lp_y = df['ne (1e18 m-3)'].values[:250]
bin_lp_x = df['R-Rsep OMP (cm).1'].values[:250]
bin_lp_x_err = df['Error (cm)'].values[:250]
bin_lp_y = df['ne (1e18 m-3).1'].values[:250]
bin_lp_y_err = df['Error'].values[:250]
bin_lp_x = bin_lp_x[~np.isnan(bin_lp_x)][2:]
bin_lp_y = bin_lp_y[~np.isnan(bin_lp_y)][2:]
bin_lp_x_err = bin_lp_x_err[~np.isnan(bin_lp_x_err)][2:]
bin_lp_y_err = bin_lp_y_err[~np.isnan(bin_lp_y_err)][2:]


fig = pp.pplot(all_lp_x, all_lp_y, alpha=0.1, ms=5, logy=True, color=20)
fig = pp.pplot(bin_lp_x, bin_lp_y, logy=True, yrange=[0.1, 10], xrange=[4, 14], fig=fig, color=20, xerr=bin_lp_x_err, yerr=bin_lp_y_err)

main_exp_a = 25.185
main_exp_b = -0.252
wind_exp_a = 23880
wind_exp_b = -0.845

main_exp_x = np.linspace(4, 11.5, 100)
wind_exp_x = np.linspace(11.787, 14, 100)
main_exp_y = main_exp_a * np.exp(main_exp_b * main_exp_x)
wind_exp_y = wind_exp_a * np.exp(wind_exp_b * wind_exp_x)

fig = pp.pplot(main_exp_x, main_exp_y, fmt='--', fig=fig, color=20, lw=3)
fig = pp.pplot(wind_exp_x, wind_exp_y, fmt='--', fig=fig, color=20, lw=3, xlabel='R-Rsep OMP (cm)', ylabel='ne (1e18 m-3)')

xl_path2 = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/07/LPData.xlsx'
df2 = pd.read_excel(xl_path2, sheet_name='Lengths', skiprows=2)

conn_x = df2['A.2'].values
conn_y = df2['A.6'].values

# Get the Axes object so we can twin it for a secondary y axis.
ax2 = fig.axes[0].twinx()
ax2.plot(conn_x, conn_y, lw=4, color=pp.tableau20[18])
ax2.set_ylim([0, 8])
ax2.tick_params(axis='both', which='both', labelsize=18, color=pp.tableau20[18], labelcolor=pp.tableau20[18])
ax2.set_ylabel('Connection Length (m)', color=pp.tableau20[18], fontsize=26)
fig.tight_layout()


# Create another plot, just with the A2 data instead of LP.
xl_path3 = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx'
df3 = pd.read_excel(xl_path3)

a2_x = df3['R-Rsep omp D (cm)'].values[10:]
a2_y = df3['W Areal Density D (1e15 W/cm2)'].values[10:]
a2_x_err = df3['R-Rsep omp Error D (cm)'].values[10:]
a2_y_err = df3['W Areal Density Error D (1e15 W/cm2)'].values[10:]

fig = pp.pplot(a2_x, a2_y, xerr=a2_x_err, yerr=a2_y_err, color=20, xrange=[4, 14], logy=True)

a2_main_exp_a = 2.0034
a2_main_exp_b = -0.201
a2_wind_exp_a = 5.4521E6
a2_wind_exp_b = -1.719

a2_main_exp_x = np.linspace(7, 9.5, 100)
a2_main_exp_y = a2_main_exp_a * np.exp(a2_main_exp_b * a2_main_exp_x)
a2_wind_exp_x = np.linspace(10, 12, 100)
a2_wind_exp_y = a2_wind_exp_a * np.exp(a2_wind_exp_b * a2_wind_exp_x)

fig = pp.pplot(a2_main_exp_x, a2_main_exp_y, color=20, fmt='--', lw=3, fig=fig)
fig = pp.pplot(a2_wind_exp_x, a2_wind_exp_y, color=20, fmt='--', lw=3, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)')

# Get the Axes object so we can twin it for a secondary y axis.
ax2 = fig.axes[0].twinx()
ax2.plot(conn_x, conn_y, lw=4, color=pp.tableau20[18])
ax2.set_ylim([0, 8])
ax2.tick_params(axis='both', which='both', labelsize=18, color=pp.tableau20[18], labelcolor=pp.tableau20[18])
ax2.set_ylabel('Connection Length (m)', color=pp.tableau20[18], fontsize=26)
fig.tight_layout()
