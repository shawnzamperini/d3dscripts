import pretty_plots as pp
import numpy as np
import pandas as pd


a8_itf_r = np.array([18.63, 18.04, 17.45, 16.87, 16.29, 15.71, 15.15, 14.58,
                     14.02, 13.46, 12.90, 12.36, 11.81, 11.26, 10.72, 10.18,
                     9.65,  9.12,  8.59,  8.06])
a8_otf_r = np.array([19.66, 19.06, 18.46, 17.87, 17.29, 16.71, 16.13, 15.55,
                     14.99, 14.42, 13.86, 13.30, 12.75, 12.20, 11.65, 11.11,
                     10.57, 10.03, 9.50,  8.97])
a8_itf_w = np.array([0.0000, 0.0000, 0.0000, 0.0000,0.0000,0.0014,0.0014,0.0000,
0.0000,0.0014,0.0000,0.0055,0.0138,0.0165,0.0193,0.0289,0.0234,0.0206,0.0303,0.0069])
a8_otf_w = np.array([0.0000,0.0000,0.0000,0.0000,0.0000,0.0014,0.0000,0.0000,0.0000,
0.0014,0.0014,0.0000,0.0028,0.0014,0.0069,0.0041,0.0124,0.0069,0.0193,0.0014])
a15_itf_r = np.array([17.75,17.17,16.59,16.02,15.45,14.89,14.33,13.77,13.22,12.67,
12.12,11.58,11.04,10.50,9.97,9.44,8.91,8.39,7.86,7.35])
a15_otf_r = np.array([16.75,16.17,15.60,15.04,14.48,13.92,13.37,12.82,12.27,
11.73,11.18,10.65,10.11,9.58,9.06,8.53,8.01,7.49,6.97,6.46])
a15_itf_w = np.array([0.0014,0.0000,0.0000,0.0000,0.0014,0.0000,0.0000,0.0000,
0.0000,0.0014,0.0000,0.0014,0.0000,0.0027,0.0082,0.0068,0.0041,0.0000,0.0000,0.0164])
a15_otf_w = np.array([0.0000,0.0041,0.0014,0.0000,0.0000,0.0041,0.0027,0.0041,
0.0000,0.0041,0.0041,0.0014,0.0123,0.0123,0.0068,0.0000,0.0000,0.0205,0.0164,0.0136])


fig = pp.pplot(a8_itf_r[:-1], a8_itf_w[:-1], color=8, label='R-ITF', fmt='-')
fig = pp.pplot(a8_otf_r[:-1], a8_otf_w[:-1], color=8, label='R-OTF', fmt='--', fig=fig)
fig = pp.pplot(a15_itf_r, a15_itf_w, color=6, label='F-ITF', fmt='-', fig=fig)
fig = pp.pplot(a15_otf_r, a15_otf_w, color=6, label='F-OTF', fmt='--', fig=fig,
               xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)', xrange=[6,18])


filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/My Slides and Sheets/2019/03/aprobe_totals.xlsx'
df = pd.read_excel(filename, sheet_name='Num. Lambdas Along', header=1)
print("Double check you are using the right columns for A8 and A15 future Shawn.")

a8_along_x = df["# of λ's.4"].values
a8_along_ratio = df['ITF/OTF.4'].values
a15_along_x = df["# of λ's.5"].values
a15_along_ratio = df['ITF/OTF.5'].values

fig = pp.pplot(a8_along_x, a8_along_ratio, fmt='-', color=8)
fig = pp.pplot(a15_along_x, a15_along_ratio, fmt='-', color=6, fig=fig,
               xlabel='# of lambdas from separatrix', ylabel='ITF/OTF Ratio')
