import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from scipy.interpolate   import griddata
from scipy.interpolate   import spline
from scipy.signal        import savgol_filter
from matplotlib          import cm


probe = 'B'
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
Z = np.clip(Z, 0, 20000)

# Centerline data. Grab it where y = 2.540275.
if probe == 'B':
    idx = np.where(y==2.540275)
elif probe == 'C':
    idx = np.where(y==2.03265)
x_center = x[idx]
z_center = z[idx]

# Edge data.
if probe == 'B':
    idx = np.where(y==0)
x_edge = x[idx]
z_edge = z[idx]

# Averaged across width of probe.
xpoints = np.unique(x).size
ypoints = np.unique(y).size
z_avgs = np.array([])
for loc in range(0, xpoints):
    avg = np.mean(z[loc*ypoints:(loc+1)*ypoints])
    z_avgs = np.append(z_avgs, avg)

# Figure of contour map.
if False:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
    cont = ax1.contourf(X/10, Y, Z, cmap='Reds')
    fig.colorbar(cont).ax.set_ylabel('LAMS Counts')
    ax1.set_xlabel('r (cm)')
    ax1.set_ylabel('z (mm)')
    fig.show()

# Figure of centerline with edge over it.
if True:
    plt.style.use('seaborn')
    fig = plt.figure()
    ax2 = fig.add_subplot(111)

    # Smoothe the data.
    x_smooth_cen = np.linspace(x_center.min(), x_center.max(), 1000)
    x_smooth_edg = np.linspace(x_edge.min(), x_edge.max(), 1000)
    z_smooth_cen = spline(x_center, z_center, x_smooth_cen)
    z_smooth_edg = spline(x_edge, z_edge, x_smooth_edg)
    z_filt_cen   = savgol_filter(z_smooth_cen, 91, 2)
    z_filt_edg   = savgol_filter(z_smooth_edg, 91, 2)
    ax2.plot(x_smooth_cen/10, z_filt_cen, label='Centerline', color='C0')
    ax2.plot(x_smooth_edg/10, z_filt_edg, label='Edge', color='C2')
    ax2.set_xlabel('Axial (cm)', fontsize=18)
    ax2.set_ylabel('LAMS Counts', fontsize=18)
    ax2.set_title('Axial Profiles of B6-ITF', fontsize=22)
    ax2.legend(frameon=True, fontsize=18)
    fig.show()

# Figure of averages.
if False:
    fig = plt.figure()
    ax3 = fig.add_subplot(111)
    ax3.plot(x_center/10, z_avgs)
    ax3.set_xlabel('r (cm)')
    ax3.set_ylabel('LAMS Counts')
    ax3.set_title('Average Profile')
    fig.show()
