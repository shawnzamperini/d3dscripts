import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from scipy.interpolate   import griddata
from matplotlib          import cm


probe = 'C'
# First load the Excel data into DataFrames.
if probe == 'B':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/BU6_Map.xlsx',
                       sheet_name = 'Sheet3',
                       usecols    = [0, 2, 6],
                       names      = ['Z', 'X', 'Y'])
elif probe == 'C':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/CU6_Map.xlsx',
                       sheet_name = 'Sheet2',
                       usecols    = [0, 2, 6],
                       names      = ['Z', 'X', 'Y'])

# Put the data in 2D formats for countour plotting.
x = df['X'].values
y = df['Y'].values
z = df['Z'].values
X, Y = np.meshgrid(x, y)
Z = griddata(points = (x, y),
             values = z,
             xi     = (X, Y),
             method = 'nearest')

# Centerline data. Grab it where y = 2.540275.
if probe == 'B':
    idx = np.where(y==2.540275)
elif probe == 'C':
    idx = np.where(y==2.03265)
x_center = x[idx]
z_center = z[idx]

# Averaged across width of probe.
xpoints = np.unique(x).size
ypoints = np.unique(y).size
z_avgs = np.array([])
for loc in range(0, xpoints):
    avg = np.mean(z[loc*ypoints:(loc+1)*ypoints])
    z_avgs = np.append(z_avgs, avg)

# Figure of contour map.
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
cont = ax1.contourf(X/10, Y, Z, cmap='Reds')
fig.colorbar(cont).ax.set_ylabel('LAMS Counts')
ax1.set_xlabel('r (cm)')
ax1.set_ylabel('z (mm)')
fig.show()

# Figure of centerline.
fig = plt.figure()
ax2 = fig.add_subplot(111)
ax2.plot(x_center/10, z_center)
ax2.set_xlabel('r (cm)')
ax2.set_ylabel('LAMS Counts')
ax2.set_title('Centerline Profile')
fig.show()

# Figure of averages.
fig = plt.figure()
ax3 = fig.add_subplot(111)
ax3.plot(x_center/10, z_avgs)
ax3.set_xlabel('r (cm)')
ax3.set_ylabel('LAMS Counts')
ax3.set_title('Average Profile')
fig.show()
