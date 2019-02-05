import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from matplotlib          import cm

max_val = 20000
filename = '/home/shawn/d3dscripts/Data/Polodial_Scans/BU06_Map_125um_spacing_New.xlsx'
df = pd.read_excel(filename, sheet_name='Sheet1')

df_tot = pd.DataFrame()
for pol_loc in np.arange(0, 5.0, 0.125):
    idx = np.where(df['Poloidal [mm]'].values == pol_loc)
    axial_locs = df['Radial [mm]'].values[idx]
    df_tot[pol_loc] = df['total'].values[idx]
df_tot.set_index(np.unique(df['Radial [mm]']), inplace=True)

df_182 = pd.DataFrame()
for pol_loc in np.arange(0, 5.0, 0.125):
    idx = np.where(df['Poloidal [mm]'].values == pol_loc)
    axial_locs = df['Radial [mm]'].values[idx]
    df_182[pol_loc] = df['EF182'].values[idx]
df_182.set_index(np.unique(df['Radial [mm]']), inplace=True)

# Get in numpy array form for easier manipulation.
arr_tot = df_tot.as_matrix()
arr_182 = df_182.as_matrix()

# Anything with less then 0.6 182EF, set to zero (for now).
mask1 = arr_182 < 0.0
arr_tot[mask1] = 0.0

# Create grid for the countour plot.
x = df_tot.index.values
y = df_tot.columns.values
Y, X = np.meshgrid(y, x)

# Replace zeros with averages of the column (i.e. the axial location).
for col in np.arange(0, x.size):
    # Get non zero values
    mask2 = arr_tot[:][col] > 0.0

    # Get the average of the non zero values.
    avg = np.mean(arr_tot[:][col][mask2])

    # Replace the zeros with the average.
    arr_tot[:][col][np.logical_not(mask2)] = avg

Z = arr_tot

if max_val != None:
    Z = np.clip(arr_tot, 0, max_val)  # To prevent errant points ruining the countour plot.

# Plot 2D contour map.
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
cont = ax1.contourf(X/10, Y, Z, cmap='Reds')
fig.colorbar(cont).ax.set_ylabel('LAMS Counts')
ax1.set_xlabel('Axial Length (cm)')
ax1.set_ylabel('Poloidal Length (mm)')
fig.show()
