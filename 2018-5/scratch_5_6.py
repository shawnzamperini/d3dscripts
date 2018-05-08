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

    # Get data into variables.
    x = df['X'].values[:2478] # ignore the last point at 60 since its different amount of measurements.
    y = df['Y'].values[:2478]
    z = df['Z'].values[:2478]

    xs = np.unique(x)
    ys = np.unique(y)[:-2]
    data = z.reshape(xs.size, ys.size)
    lams = pd.DataFrame(data, index=xs, columns=ys)
    avg_lams = pd.DataFrame()

    for i in range(0, 14):
        avg_lams['loc ' + str(i)] = lams.iloc[i*4:(i+1)*4].mean()

    # Divide axial direction (x) into 15 segments (60 locations total), so every 4 x points.
    x_avgs = np.array([])
    for i in range(0, 15):
        x_vals = np.unique(x)[4*i:4*(i+1)]
        x_avg = np.mean(x_vals)
        x_avgs = np.append(x_avgs, x_avg)

    avg_lams.columns = x_avgs[:-1] # Bc we're ignoring the last data point.

elif probe == 'C':
    df = pd.read_excel(io         = '/home/shawn/d3dscripts/Data/Polodial_Scans/CU6_Map.xlsx',
                       sheet_name = 'Sheet2',
                       usecols    = [0, 2, 6],
                       names      = ['Z', 'X', 'Y'])

    x = df['X'].values
    y = df['Y'].values
    z = df['Z'].values

    xs = np.unique(x)
    ys = np.unique(y)
    data = z.reshape(xs.size, ys.size)
    lams = pd.DataFrame(data, index=xs, columns=ys)
    avg_lams = pd.DataFrame()

    for i in range(0, 14):
        avg_lams['loc ' + str(i)] = lams.iloc[i*4:(i+1)*4].mean()

    # Divide axial direction (x) into 15 segments (60 locations total), so every 4 x points.
    x_avgs = np.array([])
    for i in range(0, 15):
        x_vals = np.unique(x)[4*i:4*(i+1)]
        x_avg = np.mean(x_vals)
        x_avgs = np.append(x_avgs, x_avg)

    avg_lams.columns = x_avgs[:-1] # Bc we're ignoring the last data point.
