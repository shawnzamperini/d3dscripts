import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from scipy.interpolate   import griddata
from matplotlib          import cm
import matplotlib as mpl

df_tot = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map_New.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:W')
df_tot = df_tot.iloc[0:202]
df_tot = df_tot.drop('Poloidal [mm]', axis=1)
max_val = 20000

df_182 = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map_New.xlsx', sheet_name='EF Map', index_col=0, usecols='A:W')
df_182 = df_182.iloc[0:202]
df_182 = df_182.drop('Poloidal [mm]', axis=1)
#max_val = 0.95

arr_182 = df_182.as_matrix()
arr_tot = df_tot.as_matrix()

mask = arr_182 < 0.6
arr_tot[mask] = 0.0

x = df_tot.index.values
y = df_tot.columns.values
Y, X = np.meshgrid(y, x)

# Replace zeros with averages of the column (i.e. the axial location).
for col in np.arange(0, x.size):
    # Get non zero values
    mask = arr_tot[:][col] > 0.0

    # Get the average of the non zero values.
    avg = np.mean(arr_tot[:][col][mask])

    # Replace the zeros with the average.
    arr_tot[:][col][np.logical_not(mask)] = avg

Z = arr_tot

if max_val != None:
    Z = np.clip(arr_tot, 0, max_val)  # To prevent errant points ruining the countour plot.

norm = mpl.colors.Normalize(vmin=0, vmax=max_val)

# Plot 2D contour map.
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
cont = ax1.contourf(X/10, Y, norm(Z), cmap='Reds')
cbar = fig.colorbar(cont, ticks=np.linspace(0,1,6))
cbar.ax.set_ylabel('LAMS Counts', size=24, weight='bold')
cbar.ax.tick_params(labelsize=20)
ax1.set_xlabel('Axial Length (cm)', size=24, weight='bold')
ax1.set_ylabel('Poloidal Length (mm)', size=24, weight='bold')
ax1.tick_params(labelsize=22)
fig.tight_layout()
fig.show()
