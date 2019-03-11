import pandas       as pd
import pretty_plots as pp


path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A31.xlsx"
df = pd.read_excel(path)

itf_x = df['R-Rsep omp D (cm)'].values[:-1]
itf_y = df['W Areal Density D (1e15 W/cm2)'].values[:-1]
otf_x = df['R-Rsep omp U (cm)'].values[:-1]
otf_y = df['W Areal Density U (1e15 W/cm2)'].values[:-1]

fig = pp.pplot(itf_x, itf_y, fmt='-', label='ITF')
fig = pp.pplot(otf_x, otf_y, fmt='-', label='OTF', fig=fig, color=8,
               xlabel='R-Rsep omp (cm)', ylabel='W Areal Density (1e15 W/cm2)')
