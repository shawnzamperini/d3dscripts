import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pretty_plots as pp
import sys

sys.path.append('/mnt/c/Users/Shawn/Documents/GitHub/d3dscripts/2019/09')
from sim_probes_v2 import plot_with_rbs

a8_rbs = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
         'Collector Probe Excel Sheets/A8.xlsx'
a15_rbs = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
         'Collector Probe Excel Sheets/A15.xlsx'

bd10_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/BD10_Map_Analysis.xlsx'
au08_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AU08_Map_Analysis.xlsx'
ad08_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AD08_Map_Analysis.xlsx'
au15_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AU15_Map_Analysis.xlsx'
ad15_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AD15_Map_Analysis.xlsx'

# Some constants.
middle = 2.0
avg = True
ignore = 15

# Run above function to get the data.
a08itf = plot_with_rbs(ad08_file, a8_rbs, 'D', 0.5E-06,   0, color=8, middle=middle, avg=avg)
a08otf = plot_with_rbs(au08_file, a8_rbs, 'U', 0.5E-06,   0, color=8, middle=middle, avg=avg)
a15itf = plot_with_rbs(au15_file, a15_rbs, 'U', 5.015E-07, 0, r_shift=0, middle=middle, avg=avg)
a15otf = plot_with_rbs(ad15_file, a15_rbs, 'D', 5.015E-07, 0, r_shift=0, middle=middle, avg=avg)

# Some errant data points.
a08itf['LAMS W'][222] = a08itf['LAMS W'][221]
a08itf['LAMS W Error'][222] = a08itf['LAMS W Error'][221]

fig, ax = plt.subplots()

# A8
ignore = 15
ax.fill_between(a08itf['LAMS Romp'][ignore:], a08itf['LAMS W'][ignore:]-a08itf['LAMS W Error'][ignore:],
                a08itf['LAMS W'][ignore:]+a08itf['LAMS W Error'][ignore:], color=pp.tableau20[18],
                alpha=0.25)
ax.plot(a08itf['LAMS Romp'][ignore:], a08itf['LAMS W'][ignore:], color=pp.tableau20[18], lw=2, label='LAMS')

# A15
ignore = 0
ax.fill_between(a15itf['LAMS Romp'][ignore:], a15itf['LAMS W'][ignore:]-a15itf['LAMS W Error'][ignore:],
                a15itf['LAMS W'][ignore:]+a15itf['LAMS W Error'][ignore:], color=pp.tableau20[12],
                alpha=0.25)
ax.plot(a15itf['LAMS Romp'][ignore:], a15itf['LAMS W'][ignore:], color=pp.tableau20[12], lw=2, label='LAMS')
ax.errorbar(a08itf['RBS Romp'][:-1], a08itf['RBS W'][:-1], yerr=0.004, ecolor='k',
            color=pp.tableau20[18], marker='.', mec='k', ms=10, linestyle='None', label='RBS')
ax.errorbar(a15itf['RBS Romp'], a15itf['RBS W'], yerr=0.004, ecolor='k',
            color=pp.tableau20[12], marker='.', mec='k', ms=10, linestyle='None', label='RBS')
ax.set_ylim([0, None])
ax.set_xlim([6, 15])
ax.legend(fontsize=16)
ax.set_xlabel('R-Rsep OMP (cm)', fontsize=16)
ax.set_ylabel('W Areal Density (1e15 W/m2)', fontsize=16)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
