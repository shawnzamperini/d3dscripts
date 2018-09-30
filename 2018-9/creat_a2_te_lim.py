import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

filename = '/mnt/c/Users/Shawn/Documents/d3d_work/A2_te.xlsx'
df = pd.read_excel(filename, usecols=[0,1,2,3])

def exp_fit(x, a, b):
    return a * np.exp (-x * b)

popt_te, pcov_te = curve_fit(exp_fit, df['R in LIM Coordinates (m)'], df['Fit Te (eV)'])
popt_ne, pcov_ne = curve_fit(exp_fit, df.dropna()['R in LIM Coordinates (m)'], df.dropna()['ne'])
lim_rs = np.arange(-0.15, 0.062, 0.002)
lim_te = exp_fit(lim_rs, *popt_te)
lim_ne = exp_fit(lim_rs, *popt_ne)

#for r, te in list(zip(lim_rs, lim_te)):
#    print("{:.3f}, {:.3f}".format(r, te))

for r, ne in list(zip(lim_rs, lim_ne)):
    print("{:.3f} {:.3e}".format(r, ne*10**18))
