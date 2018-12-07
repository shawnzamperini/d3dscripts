import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib          import cm
import matplotlib        as mpl
import pretty_plots as pp


filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/CU07_Map.xlsx'

df = pd.read_excel(filename, sheet_name='Sheet2')[:-1].set_index('Radial [mm]')

X, Y = np.meshgrid(df.index.values, df.columns.values[:-1][::-1])
#Z = df.values.T

max_val = 100000
Z = np.clip(df.values[:,:-1], 0, max_val).T
norm = mpl.colors.Normalize(vmin=0, vmax=max_val)

fig = plt.figure(figsize=(10, 7.5))
ax1 = fig.add_subplot(111)
fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
cont = ax1.contourf(X/10, Y, norm(Z), cmap='Reds')
cbar = fig.colorbar(cont, ticks=np.linspace(0,1,6))
cbar.ax.set_ylabel('LAMS Counts', size=24, weight='bold')
cbar.ax.tick_params(labelsize=20)
ax1.set_xlabel('Axial Length (cm)', size=24, weight='bold')
ax1.set_ylabel('Z Location (mm)', size=24, weight='bold')
ax1.tick_params(labelsize=22)

# Throw each measurement location on.
ax1.plot(X/10, Y, 'k.', ms=3)

fig.tight_layout()
fig.show()

# Make plot of the average poloidal profile.
avg_pol_counts = pd.read_excel(filename, sheet_name='Sheet2').set_index('Radial [mm]').loc['Average']
x = avg_pol_counts.index.values[:-1][::-1]
y = avg_pol_counts.values[:-1] / np.max(avg_pol_counts.values[:-1])
pp.pplot(x, y, fmt='-', xlabel='Z Location (mm)', ylabel='Axial Averaged Counts', weight='bold')
