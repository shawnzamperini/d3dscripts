import pretty_plots as pp
import pandas as pd


path = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/My Slides and Sheets/2019/03/floor_frac.xlsx'
df = pd.read_excel(path, sheet_name='Sheet3', skiprows=16, usecols=[0, 1], names=['Probe', 'Floor Frac.'])

pnames = df['Probe'].values
fracs  = df['Floor Frac.'].values

fig = pp.pbar(fracs[:-4], y2=fracs[-4:], ylabel='Floor Fraction', label='Forward', label2='Reverse',
        bar_names=pnames[:-4], bar_names2=[pnames[-4:]], tick_rot=-90)
