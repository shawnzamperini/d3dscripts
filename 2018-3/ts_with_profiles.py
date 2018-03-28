import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import get_ts as ts
import MDSplus as mds
import scipy.interpolate as interpolate
from scipy.optimize import curve_fit

# Load TS data for shot 167405.
ts_df = pd.read_excel('/home/shawn/d3dscripts/Data/ts167405.xlsx')
# Only want SOL data (psin > 1.0).
sol_df = ts_df.iloc[:16].copy()

# Already did this and coped the data to the Excel.
if False:
    # Go from TS psin -> R-Rsep omp. First MDSplus connection.
    conn     = mds.Connection("localhost")
    # Then load gfile from where plasma is stable.
    gfile    = ts.load_gfile_mds(167405, 3000, connection=conn)
    # Get R, Z grid.
    Rs, Zs   = np.meshgrid(gfile['R'], gfile['Z'])
    # Get R and Z of magentic axis.
    Z_axis   = gfile['ZmAxis']
    R_axis   = gfile['RmAxis']
    # Get R's and Z's of separatrix.
    Zes      = np.copy(gfile['lcfs'][:, 1][13:-17])
    Res      = np.copy(gfile['lcfs'][:, 0][13:-17])
    # Just get right half of the coordinates.
    Rs_trunc = Rs > R_axis
    # Interpolations for Rsep(Z), Psin(R, Z), and Romp(psin, Z).
    f_Rs     = interpolate.interp1d(Zes, Res, assume_sorted=False)
    f_psin   = interpolate.Rbf(Rs, Zs, gfile["psiRZn"])
    f_Romp   = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
    # R of separatrix at omp.
    rsep_omp = f_Rs(Z_axis)

    # Create R-Rsep omp for TS psins.
    omp_arr = np.array([])
    for psin in sol_df['Psin']:
        omp_arr = np.append(omp_arr, f_Romp(psin, Z_axis) - rsep_omp)
    sol_df['R-Rsep omp'] = omp_arr

# Now get the probe data for B10/C10.
b10_df = pd.read_excel('/home/shawn/d3dscripts/Data/B10.xlsx', header=1)
c10_df = pd.read_excel('/home/shawn/d3dscripts/Data/C10.xlsx', header=1)

def exp_fit(x, a, b, c):
    return a * np.exp(-x * b) + c

tofitx = np.array(sol_df['R-Rsep omp'] * 100, dtype=np.float64)
tofity = np.array(sol_df['Flux']*10**-21, dtype=np.float64)

guess = (200, 60, 0)
popt, pcov = curve_fit(exp_fit, tofitx, tofity, p0=guess, maxfev=50000)
x_range = np.linspace(0, 20, 100)
y_fitted = exp_fit(x_range, *popt)

fig1, ax1 = plt.subplots()
ax1.plot(sol_df['R-Rsep omp'] * 100.0, sol_df['Flux']*10**-21, 'r')
ax1.plot(x_range, y_fitted, 'r--', label='Exp. Fit')
ax1.set_xlabel('R-Rsep omp')
ax1.set_ylabel('D Flux (10**21 m-2 s-1)', color='r')
ax1.tick_params('y', colors='red')
ax2 = ax1.twinx()
ax2.plot(b10_df['rminrsep_omp_U'], b10_df['w_areal_U'], 'bv', label='B10')
ax2.plot(c10_df['rminrsep_omp_U'], c10_df['w_areal_U'], 'b1', label='C10')
ax2.tick_params('y', colors='blue')
plt.legend()
plt.show()
