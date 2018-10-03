import sys
sys.path.append('/home/shawn/d3dscripts/Data/Polodial_Scans/')

import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from scipy.interpolate   import griddata
from matplotlib          import cm
import matplotlib        as mpl


probe = 'BD06'

if probe == 'BD06':
    # Load data into DataFrame.
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map.xlsx', sheet_name='Sheet1', index_col=0)

    # Drop the table at the end.
    df = df.iloc[0:202]

    # Drop the column of all zeros.
    df = df.drop('Poloidal [mm]', axis=1)

    # Cap max size at 20,000.
    max_val = 20000

elif probe == 'CU06':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/CU6_Map.xlsx',
                       sheet_name = 'Sheet2',
                       usecols    = [0, 2, 6],
                       names      = ['Z', 'X', 'Y'])
    x = df['X'].values
    y = df['Y'].values
    z = df['Z'].values
    df = pd.DataFrame(z.reshape(np.unique(x).size, np.unique(y).size))
    df.index   = np.unique(x)
    df.columns = np.unique(y)
    max_val = 20000
elif probe == 'CD06':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/CD06_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:M')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 18000
elif probe == 'AU03':
    df = pd.read_excel('AU03_Map.xlsx', sheet_name='Sheet1', index_col=0)
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 150000
elif probe == 'AD03':
    df = pd.read_excel('AD03_Map.xlsx', sheet_name='Sheet1', index_col=0)
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 800000
elif probe == 'AD04':
    df = pd.read_excel('AD04_Map.xlsx', sheet_name='Sheet1', index_col=0)
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 100000
elif probe == 'AU04':
    df = pd.read_excel('AU04_Map.xlsx', sheet_name='Sheet1', index_col=0)
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 50000
elif probe == 'AD09':
    df = pd.read_excel('AD09_Map.xlsx', sheet_name='Sheet1', index_col=0)
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 25000
elif probe == 'AU09':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/AU9_Map.xlsx',
                       sheet_name = 'Sheet2',
                       usecols    = [0, 2, 6],
                       names      = ['Z', 'X', 'Y'])
    x = df['X'].values[19:]
    y = df['Y'].values[19:]
    z = df['Z'].values[19:]
    df = pd.DataFrame(z.reshape(np.unique(x).size, np.unique(y).size))
    df.index   = np.unique(x)
    df.columns = np.unique(y)
    max_val = None
elif probe == 'CU04':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/CU04_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:N')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 150000
elif probe == 'BD06_New_Total':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map_New.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 20000
elif probe == 'BD06_New_182': # Note this is for the EF map, not counts.
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map_New.xlsx', sheet_name='EF Map', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 0.9
elif probe == 'BU04':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BU04_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 400000
elif probe == 'BD04':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/BD04_Map.xlsx',
                       sheet_name = 'Sheet1',
                       usecols    = [0, 1, 2],
                       names      = ['X', 'Y', 'Z'])
    y = df['X'].values
    x = df['Y'].values
    z = df['Z'].values
    df = pd.DataFrame(np.transpose(z.reshape(np.unique(x).size, np.unique(y).size))) # Note have to transpose and swap x, y bc the data is shaped differently.
    df.index   = np.unique(y)
    df.columns = np.unique(x)
    max_val = 500000
elif probe == 'AU06':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/AU06_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:T')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 50000
elif probe == 'AD06':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/AD06_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:T')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 165000
elif probe == 'BD03':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD03_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 2400000
elif probe == 'BU03':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BU03_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 2000000
elif probe == 'CU03':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/CU03_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:M')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 250000
elif probe == 'CD03':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/CD03_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:M')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 600000
elif probe == 'AD05':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/AD05_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:T')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 500000
elif probe == 'AU05':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/AU05_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:T')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 100000
elif probe == 'CD04':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/CD04_Map.xlsx',
                       sheet_name = 'Sheet1',
                       usecols    = [0, 1, 2],
                       names      = ['X', 'Y', 'Z'])
    y = df['X'].values
    x = df['Y'].values
    z = df['Z'].values
    df = pd.DataFrame(np.transpose(z.reshape(np.unique(x).size, np.unique(y).size))) # Note have to transpose and swap x, y bc the data is shaped differently.
    df.index   = np.unique(y)
    df.columns = np.unique(x)
    max_val = 300000


# Grid the data for 2D plotting. Can use max_val as either max counts or max enrichment fraction!
x    = df.index.values
y    = df.columns.values
Y, X = np.meshgrid(y, x)
if max_val != None:
    Z = np.clip(df.as_matrix(), 0, max_val)  # To prevent errant points ruining the countour plot.
else:
    Z = df.as_matrix()

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
