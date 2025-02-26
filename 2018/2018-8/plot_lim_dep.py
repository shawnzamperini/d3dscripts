import matplotlib.pyplot  as plt
import pandas             as pd
import numpy              as np
from mpl_toolkits.mplot3d import Axes3D
import pretty_plots as pp


# Location of Excel file with "TOTAL DEPOSITION" table copy/pasted into it.
filename = '/mnt/c/Users/Shawn/Documents/d3d_work/3DLIM Runs/test22b_totdep.xlsx'

# Load the data into a Dataframe, drop the last row that has junk in it, and
# then rename the columns that is a string '0.1' to a float 0.0 (not sure why
# this happens).
print("Loading file...")
df = pd.read_excel(filename, index_col=0)
df.drop(df.columns[-1], axis=1, inplace=True)
df.rename({'0.1':0.0, '0.2':0.0}, axis=1, inplace=True)

# Create the meshgrid for use in the 3D plots (X = poloidal, Y = radial?
# Z = counts). Z is multiplied by -1 because they're counts of erosion, so
# deposition shows up as a negative number.
print("Creating mesh...")
X, Y = np.meshgrid(np.float64(df.columns), np.float64(df.index))
# Only plot one side.
X, Y = np.meshgrid(np.float64(df.columns), np.float64(df.index[df.index>0]))
Z = df.iloc[df.index>0].values * -1

# Plotting commands.
def plot_3d(angle1, angle2):
    font = {'fontsize':18}
    fig = plt.figure(figsize=((15,8)))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='Reds')
    ax.set_xlabel('\nPoloidal (m)', font)
    ax.set_ylabel('\nDistance along probe (m)', font)
    ax.set_zlabel('\nDeposition counts', font)
    ax.view_init(angle1, angle2)
    #fig.tight_layout()
    fig.show()

# Get each side of probe, and a centerline (P=0).
itf_x = df.iloc[df.index > 0.0][0].index.values
itf_y = df.iloc[df.index > 0.0][0].values * -1
otf_x = df.iloc[df.index < 0.0][0].index.values * -1
otf_y = df.iloc[df.index < 0.0][0].values * -1

fig = pp.pplot(itf_x, itf_y, fmt='-', label='ITF')
fig = pp.pplot(otf_x, otf_y, fmt='-', fig=fig, color=8, label='OTF',
               xlabel='Distance along probe (m)',
               ylabel='Deposition (arbitrary units)')

print("Creating plots...")
plot_3d(35, -85)
plot_3d(35, -32)
print("Done.")
