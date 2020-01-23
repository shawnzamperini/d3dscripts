import pandas as pd
import numpy  as np
import pretty_plots as pp


xl_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/11/match_3dlim_real.xlsx"

df = pd.read_excel(xl_path, sheet_name='Python Sheet')

lw=3
ms=14
mew=3

fig = pp.pplot(df['ITF X'], df['ITF Y (50)'], fmt='-', lw=lw, logy=True)
fig = pp.pplot(df['ITF X'], df['ITF Y (5)'],  fmt='-', lw=lw*2, fig=fig)
fig = pp.pplot(df['OTF X'], df['OTF Y'], color=8, fmt='-', lw=lw, fig=fig)
fig = pp.pplot(df['ITF RBS X'], df['ITF RBS Y'], ms=ms, markeredgewidth=mew, logy=True, fig=fig)
fig = pp.pplot(df['OTF RBS X'], df['OTF RBS Y'], ms=ms, color=8, xrange=[7,14], markeredgewidth=mew, xlabel='R-Rsep OMP (cm)', ylabel='Deposition (Normalized)', fig=fig)
fig.axes[0].grid(True, which='major')
