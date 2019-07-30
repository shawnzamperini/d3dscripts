import pretty_plots as pp
import pandas       as pd
import numpy        as np
import matplotlib.pyplot as plt


# Which plots to plot.
plots = [3]

# Some constants for plotting.
fontsize=26

xl_path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/06/vti_and_vi.xlsx'
df1 = pd.read_excel(xl_path, sheet_name='Use OEDGE Instead', skiprows=[0,1,2], usecols=np.arange(0,10))
df2 = pd.read_excel(xl_path, sheet_name='OEDGE Mach', skiprows=1, usecols=[6,7,8,9,10])

s = df1['s (m)'].values

if 0 in plots:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(s, df1['Ti (eV)'], '-', color=pp.tableau20[6], lw=3)
    ax1.set_xlabel('s (m)', fontsize=fontsize)
    ax1.set_ylabel('Ti (eV)', color=pp.tableau20[6], fontsize=fontsize)
    ax1.tick_params(axis='y', labelcolor=pp.tableau20[6])
    ax2 = ax1.twinx()
    ax2.plot(s, df1['dTi/ds'], '-', color=pp.tableau20[8], lw=3)
    ax2.set_ylabel('dTi/ds (eV/m)', color=pp.tableau20[8], fontsize=fontsize)
    ax2.tick_params(axis='y', labelcolor=pp.tableau20[8])
    fig.tight_layout()
    fig.show()

if 1 in plots:
    fig = pp.pplot(s, df1['Ti (eV)'], fmt='-', xlabel='s (m)', ylabel='Ti (eV)')
    fig = pp.pplot(s, df1['dTi/ds'], fmt='-', xlabel='s (m)', ylabel='dTi/ds (eV/m)')

if 2 in plots:
    fig = pp.pplot(s, df2['M - R.2'], color=8, fmt='-', label='Reverse')
    fig = pp.pplot(s, df2['M - F.2'], color=6, fmt='-', label='Forward', fig=fig)

    s_look = [2.70802535, 27.0802535, 54.160507, 81.2407605, 97.4889126]
    r_look = [-0.5, -0.2, 0, 0.8, 0.5]
    f_look = [-0.5, 0.3, 0.5, 0.5, 0.3]

    fig = pp.pplot(s_look, r_look, color=8, ms=12, fig=fig)
    fig = pp.pplot(s_look, f_look, color=6, ms=12, fig=fig, xlabel='s (m)', ylabel='Mach Number')

if 3 in plots:
    fig = pp.pplot(s, df1['vTi (m/s)'], color=14, fmt='--')
    fig = pp.pplot(s, df1['vi - R'], color=8, fmt=':', label='vi', fig=fig)
    fig = pp.pplot(s, df1['vi - F'], color=6, fmt=':', label='vi', fig=fig)
    fig = pp.pplot(s, df1['vz - R'], color=8, fmt='-', label='vz', fig=fig)
    fig = pp.pplot(s, df1['vz - F'], color=6, fmt='-', label='vz', fig=fig)
