import pandas as pd
import numpy as np
import pretty_plots as pp


filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/My Slides and Sheets/2019/03/aprobe_totals.xlsx'
df = pd.read_excel(filename, sheet_name='A Probes', usecols=np.arange(0,7))[0:19]

x = df['# of lambdas'].values
y = df['ITF/OTF Total'].values
xr = np.array([]); yr = np.array([])
xf = np.array([]); yf = np.array([])
xr = np.concatenate([xr, x[:6], x[14:16]])
yr = np.concatenate([yr, y[:6], y[14:16]])
xf = np.concatenate([x[5:14], x[16:]])
yf = np.concatenate([y[5:14], y[16:]])


fig = pp.pplot(xr, yr, color=8, label='Reverse')
fig = pp.pplot(xf, yf, xlabel="# of " + r'$\mathrm{\lambda_{ne}}$' + "'s from separatrix" ,
               ylabel='Total ITF/OTF', label='Forward', fig=fig)
