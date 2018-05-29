import sys
sys.path.append('/home/shawn/d3dscripts/Data/Polodial_Scans')

import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
from scipy.interpolate   import griddata
from matplotlib          import cm


probe = 'BD06_New_182'

if probe == 'BD06':
    # Load data into DataFrame.
    df = pd.read_excel('BD06_Map.xlsx', sheet_name='Sheet1', index_col=0)

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
    df = pd.read_excel('CD06_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:M')
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
    pass
elif probe == 'CU04':
    df = pd.read_excel('CU04_Map.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:N')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 150000
elif probe == 'BD06_New_Total':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map_New.xlsx', sheet_name='Sheet1', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 20000
elif probe == 'BD06_New_182':
    df = pd.read_excel('/home/shawn/d3dscripts/Data/Polodial_Scans/BD06_Map_New.xlsx', sheet_name='EF Map', index_col=0, usecols='A:W')
    df = df.iloc[0:202]
    df = df.drop('Poloidal [mm]', axis=1)
    max_val = 0.95

# Grid the data for 2D plotting.
x    = df.index.values
y    = df.columns.values
Y, X = np.meshgrid(y, x)
if max_val != None:
    Z = np.clip(df.as_matrix(), 0, max_val)  # To prevent errant points ruining the countour plot.
else:
    Z = df.as_matrix()


# Plot 2D contour map.
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.gca().patch.set_color(cm.Reds(100))  # Fill in blank spots with low value color.
cont = ax1.contourf(X/10, Y, Z, cmap='Reds')
fig.colorbar(cont).ax.set_ylabel('LAMS Counts')
ax1.set_xlabel('Axial Length (cm)')
ax1.set_ylabel('Poloidal Length (mm)')
fig.show()
