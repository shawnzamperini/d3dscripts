import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from matplotlib          import cm

max_val = 300000
filename = '/home/shawn/d3dscripts/Data/Polodial_Scans/CD04_Map.xlsx'
df = pd.read_excel(filename, sheet_name='Sheet1')

new_df = pd.DataFrame()
for pol_loc in np.arange(0, 2.51, 0.25):
    idx = np.where(df['Poloidal [mm]'].values == pol_loc)
    axial_locs = df['Radial [mm]'].values[idx]
    new_df[pol_loc] = df['total'].values[idx]

new_df.set_index(np.unique(df['Radial [mm]']), inplace=True)

x = new_df.index.values
y = new_df.columns.values
Y, X = np.meshgrid(y, x)
if max_val != None:
    Z = np.clip(new_df.as_matrix(), 0, max_val)  # To prevent errant points ruining the countour plot.
else:
    Z = new_df.as_matrix()

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
cont = ax1.contourf(X/10, Y, Z, cmap='Reds')
fig.colorbar(cont).ax.set_ylabel('LAMS Counts')
ax1.set_xlabel('Axial Length (cm)')
ax1.set_ylabel('Poloidal Length (mm)')
fig.show()
