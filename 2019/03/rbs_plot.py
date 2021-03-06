import pandas       as pd
import pretty_plots as pp


path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
df = pd.read_excel(path)

itf_x = df['R-Rsep omp D (cm)'].values[:-3]
itf_y = df['W Areal Density D (1e15 W/cm2)'].values[:-3]
otf_x = df['R-Rsep omp U (cm)'].values[:-3]
otf_y = df['W Areal Density U (1e15 W/cm2)'].values[:-3]

fig = pp.pplot(itf_x, itf_y, fmt='-', label='ITF', lw=8, color=8)
fig = pp.pplot(otf_x, otf_y, fmt='-', label='OTF', fig=fig, color=8, lw=3,
               xlabel='R-Rsep omp (cm)', ylabel='W Areal Density (1e15 W/cm2)')
