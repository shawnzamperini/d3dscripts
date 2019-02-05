import pretty_plots as pp
import pandas as pd


leak_file = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2019-01/leakage.xlsx'
prof_file = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2019-01/plasma_profiles_from_ts.xlsx'
prob_file = '/home/shawn/Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx'

leak_df = pd.read_excel(leak_file, sheet_name='CER')
prof_df = pd.read_excel(prof_file, skiprows=[0])
prob_df = pd.read_excel(prob_file)

x1 = leak_df['Toroidal Rotation (OMP=0) (km/s)'].values
y1 = leak_df['ITF/OTF Max W Ratio'].values
x1_err = leak_df['Error.4'].values
y1_err = leak_df['Error.3'].values

"""
fig = pp.pplot(x1, y1, xerr=x1_err, yerr=y1_err,
               xlabel='Toroidal Rotation at Separatrix (km/s)',
               ylabel='Max ITF/OTF W Ratio')
"""

x2 = prof_df['Density Decay (cm)'].values
y2 = prof_df['ITF/OTF Max W Ratio'].values
x2_err = x2 * 0.1
y2_err = prof_df['Error'].values

"""
fig = pp.pplot(x2, y2, xerr=x2_err, yerr=y2_err,
               xlabel='Density Falloff (m)',
               ylabel='ITF/OTF Max W Ratio')
"""

x3 = prof_df['Density at Probe Tip'].values
x3_err = x3 * 0.1

"""
fig = pp.pplot(x3, y2, xerr=x3_err, yerr=y2_err,
               xlabel='Density at Probe Tip (m-3)',
               ylabel='ITF/OTF Max W Ratio')
"""

x4o = prob_df['R-Rsep omp U (cm)'].values
y4o = prob_df['W Areal Density U (1e15 W/cm2)'].values
x4o_err = prob_df['R-Rsep omp Error U (cm)'].values
y4o_err = prob_df['W Areal Density Error U (1e15 W/cm2)'].values

x4i = prob_df['R-Rsep omp D (cm)'].values
y4i = prob_df['W Areal Density D (1e15 W/cm2)'].values
x4i_err = prob_df['R-Rsep omp Error D (cm)'].values
y4i_err = prob_df['W Areal Density Error D (1e15 W/cm2)'].values

fig = pp.pplot(x4i, y4i, xerr=x4i_err, yerr=y4i_err, label='ITF')
fig = pp.pplot(x4o, y4o, xerr=x4o_err, yerr=y4o_err, fig=fig, label='OTF',
               xlabel='R-Rsep at OMP (cm)',
               ylabel=r'W Areal Density (1e15 W/cm$\mathrm{^2}$)',
               color=8)
