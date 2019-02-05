import pandas as pd
import pretty_plots as pp


filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
           'My Slides and Sheets/2018-10/avg_lams_profs.xlsx'

df = pd.read_excel(filename, sheet_name='Sheet3', skiprows=[0,1,2], usecols=[2,3]).dropna()
fig = pp.pplot(df['Bottom to Top (mm)'], df['Norm. Counts'], fmt='-', xlabel='', ylabel='', color=0, lw=10)
