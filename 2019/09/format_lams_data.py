import numpy as np
import pandas as pd
import pretty_plots as pp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

"""
Will just reuse the variable names from sim_probes_v2 because I'm a lazy
piece of shit.
"""

# Some constants.
middle = 2.0
avg = True
ignore = 0
lamb_ne_up = 7.81


plt.rcParams['font.family'] = 'DejaVu Sans'

# Paths to the Excel files for each side.
a08itf_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AU34_Map_Analysis.xlsx"
a08otf_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/AD34_Map_Analysis.xlsx"

# Load the RBS files.
a08rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A34.xlsx"

def plot_with_rbs(lams_path, rbs_path, u_or_d, cal_slope, cal_intercept, middle=2.0, r_shift=0, color=6, avg=False):
    """
    u_or_d: Either 'U' or 'D'.
    """

    # Load DataFrames.
    lams_df = pd.read_excel(lams_path, sheet_name='MapData')
    rbs_df  = pd.read_excel(rbs_path)

    # Get relevant RBS data.
    dist = rbs_df['Distance from Tip ' + u_or_d + ' (cm)'].values
    omp  = rbs_df['R-Rsep omp ' + u_or_d + ' (cm)'].values
    w_rbs    = rbs_df['W Areal Density ' + u_or_d + ' (1e15 W/cm2)'].values
    r_mach = rbs_df['R ' + u_or_d + ' (cm)']

    # Do a fit to get R-Rsep OMP for the LAMS data.
    p = np.polyfit(dist, omp, 1)
    lams_df['R-Rsep omp (cm)'] = p[1] + p[0] * lams_df['Axial Location [mm]'] / 10.0
    p = np.polyfit(dist, r_mach, 1)
    lams_df['R (cm)'] = p[1] + p[0] * lams_df['Axial Location [mm]'] / 10.0

    # Create a 2D DataFrame of just the Total W measurements (for RBS comparison).
    omp_locs = lams_df['R-Rsep omp (cm)'].unique()
    mach_locs = lams_df['R (cm)'].unique()
    z_locs   = lams_df['z Location [mm]'].unique()
    total_df = pd.DataFrame(); total_df2 = pd.DataFrame()  # Another for the machine R process.
    for z in z_locs:
        total_df[z]  = lams_df[lams_df['z Location [mm]'] == z]['Total W'].values
        total_df2[z] = lams_df[lams_df['z Location [mm]'] == z]['Total W'].values

    # Set index as the omp locations.
    total_df.set_index(omp_locs, inplace=True)

    # Set index of other df as the machine R locations.
    total_df2.set_index(mach_locs, inplace=True)

    if avg:
        r = total_df[middle].index.values
        r_mach_lams = total_df2[middle].index.values
        w = total_df.mean(axis=1).values
        w_err = total_df.std(axis=1).values
        w_err = w_err * cal_slope + cal_intercept

    else:
        # Get the centerline data near where RBS was taken.
        r = total_df[middle].index.values
        w = total_df[middle].values

    # Use calibration to convert from counts to areal density.
    w = w * cal_slope + cal_intercept

    # Shift LAMS R data if asked.
    r = r + r_shift

    # Plot it.
    fig = pp.pplot(r, w, fmt='-', label='LAMS', color=color)
    fig = pp.pplot(omp, w_rbs, fmt='.', ms=15, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)', label='RBS', color=color)

    # Plot the machine R plot as well.
    fig = pp.pplot(r_mach_lams, w, fmt='-', label='LAMS', color=color)
    fig = pp.pplot(r_mach, w_rbs, fmt='.', ms=15, fig=fig, xlabel='R (cm)', ylabel='W Areal Density (1e15 cm-2)', label='RBS', color=color)

    if avg:
        return {'LAMS Romp':r, 'LAMS W':w, 'LAMS W Error':w_err}
    else:
        return {'LAMS Romp':r, 'LAMS W':w}


# Run above function to get the data.
a08itf = plot_with_rbs(a08itf_path, a08rbs_path, 'U', 0.5E-06,   0, color=8, middle=middle, avg=avg)
a08otf = plot_with_rbs(a08otf_path, a08rbs_path, 'D', 0.5E-06,   0, color=8, middle=middle, avg=avg)

# Unfavorable plot.
x = a08itf['LAMS Romp'][ignore:]; y = a08itf['LAMS W'][ignore:]

# Fix errant data points.
#y[207] = y[206]
#y[301] = y[300]
#y[307] = y[306]
#y[313] = y[312]
#y[317] = y[316]

# Save for later function.
itfx_unf = x; itfy_unf = y

# Plot the ITF line.
fig = pp.pplot(x, y, fmt='-', color=12, lw=3, label='ITF')

# Add error bands if we did the average option.
if avg:
    yerr = a08itf['LAMS W Error'][ignore:]
#    yerr[207] = yerr[206]
#    yerr[301] = yerr[300]
#    yerr[307] = yerr[306]
#    yerr[313] = yerr[312]
#    yerr[317] = yerr[316]
    itfyerr_unf = yerr
    plt.fill_between(x, y - yerr, y + yerr, color=pp.tableau20[12], alpha=0.5)

# Plot the OTF line.
x = a08otf['LAMS Romp'][ignore:]; y = a08otf['LAMS W'][ignore:]

# Fix errant data point.
#y[189] = y[187]
#y[188] = y[187]

# Save for later function.
otfx_unf = x; otfy_unf = y

fig = pp.pplot(x, y, fmt='-', color=18, lw=3, fig=fig, xlabel='R-Rsep OMP (cm)', ylabel='W Areal Density (1e15 cm-2)', label='OTF', weight='bold')
plt.yticks(np.arange(0.0, 0.05, 0.01))

# Add error bands if we did the average option.
if avg:
    yerr = a08otf['LAMS W Error'][ignore:]
    yerr[189] = yerr[187]
    yerr[188] = yerr[187]
    otfyerr_unf = yerr
    plt.fill_between(x, y - yerr, y + yerr, color=pp.tableau20[18], alpha=0.5)

# Set y limits.
plt.ylim([0,None])
plt.xlabel('R-Rsep OMP (cm)')
plt.ylabel('W Areal Density (1e15 cm-2)')
plt.tight_layout()

# First create interpolation functions for the data.
f_itf_unf = interp1d(itfx_unf, itfy_unf)
f_otf_unf = interp1d(otfx_unf, otfy_unf)

# Interpolation functions for the errors as well.
f_itferr_unf = interp1d(itfx_unf, itfyerr_unf)
f_otferr_unf = interp1d(otfx_unf, otfyerr_unf)

# Find the common x range among the face (can't do ITF/OTF if there's no data!).
com_x_unf = np.linspace(np.max((itfx_unf.min(), otfx_unf.min())), np.min((itfx_unf.max(), otfx_unf.max())), 300)

# Calculate the ITF/OTF ratio over the common x range.
itf_otf_unf = f_itf_unf(com_x_unf) / f_otf_unf(com_x_unf)

# Propogate the errors.
err_unf = np.sqrt(np.power(f_itferr_unf(com_x_unf) / f_itf_unf(com_x_unf), 2) + np.power(f_otferr_unf(com_x_unf) / f_otf_unf(com_x_unf), 2))

# Plot them together.
fig = pp.pplot(com_x_unf, itf_otf_unf, fmt='-', color=8, xlabel='R-Rsep OMP (cm)', ylabel='ITF/OTF', weight='bold')

# Add error bands.
plt.fill_between(com_x_unf, itf_otf_unf - err_unf, itf_otf_unf + err_unf, color=pp.tableau20[8], alpha=0.5)

# Divided by lambda_ne plots.
fig = pp.pplot(com_x_unf / lamb_ne_up, itf_otf_unf, fmt='-', color=8, xlabel=r'R-Rsep OMP(cm) / $\mathrm{\lambda_{ne}}$', ylabel='ITF/OTF', weight='bold')

# Add error bands.
plt.fill_between(com_x_unf / lamb_ne_up, itf_otf_unf - err_unf, itf_otf_unf + err_unf, color=pp.tableau20[8], alpha=0.5)
