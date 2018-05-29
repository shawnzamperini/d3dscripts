# Add the path to where the Excel sheet is.
import sys
sys.path.append('/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets')

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize      import curve_fit


filename = '/home/shawn/Drive/School/Tennessee/Research/My Slides and Sheets/2018 - 5/max_w_vs_flux.xlsx'
df = pd.read_excel(filename, sheet_name='Ratios', usecols='A:L')

ratios      = df['Ratio'][0:11].values
ratios_err   = df['Ratio Error'][0:11].values
falloff     = df['Density Fall Off Length (m)'][0:11].values
falloff_err = df['Fall Off Error (m)'][0:11].values

# Plot the assymetry trend with the density fall off.
if True:
    font = {'fontsize' : 24,
            'weight'   : 'bold'}
    plt.style.use('seaborn')
    # Make the font larger and bold.
    #font = {'weight' : 'bold'}
    #plt.rc('font', **font)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.errorbar(falloff, ratios, ratios_err, falloff_err, 'k.', ms=20, markeredgewidth=1, capsize=5)
    ax1.set_xlabel('Density Falloff (cm)',  font)
    ax1.set_ylabel('ITF/OTF Max W Content', font)
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.tick_params(labelsize=14)
    fig.tight_layout()
    fig.show()

# Plot the asymmetry in falloff lengths.
if True:
    probe_file = '/home/shawn/d3dscripts/Data/LModeProbes.xlsx'
    aprobe_df = pd.read_excel(probe_file, sheet_name = 'A2', usecols="A:N")[:19]
    xd = aprobe_df['rminrsep_D']
    yd = aprobe_df['w_areal_D']
    xu = aprobe_df['rminrsep_U']
    yu = aprobe_df['w_areal_U']

    def exp_fit(x, a, b):
        return a * np.exp(-b * x)

    poptd, pcovd = curve_fit(exp_fit, xd[:-3], yd[:-3])
    poptu, pcovu = curve_fit(exp_fit, xu[:-3], yu[:-3])
    x_fitd = np.linspace(xd.min(), xd.max())
    x_fitu = np.linspace(xu.min(), xu.max())
    y_fitd = exp_fit(x_fitd, *poptd)
    y_fitu = exp_fit(x_fitu, *poptu)

    font = {'fontsize' : 24,
            'weight'   : 'bold'}
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(xd,     yd,     'r.', ms=20, label='A2-ITF', alpha=0.8)
    ax1.plot(x_fitd, y_fitd, 'r--', lw=3, alpha=0.8)
    ax1.plot(xu,     yu,     'b.', ms=20, label='A2-OTF', alpha=0.8)
    ax1.plot(x_fitu, y_fitu, 'b--', lw=3, alpha=0.8)
    ax1.set_xlabel(r'$\mathrm{\mathbf{R\ -\ R_{sep}\ omp\ (cm)}}$', font)
    ax1.set_ylabel(r'$\mathrm{\mathbf{W\ Density\ (10^{18}\ cm^{-2})}}$', font)
    ax1.set_ylim(top=0.6)
    ax1.tick_params(axis='both', labelsize=16)
    ax1.legend(fontsize=font['fontsize'], frameon=True, edgecolor='k')
    fig.tight_layout()
    fig.show()

# Bar graph of lambdas.
if True:
    #probes      = ['A2',  'A17',  'A18',  'A19',  'A28',  'A31',  'A32',  'A33',  'A34',  'A35']
    #lambdas_otf = [1.2547, 4.1982, 3.2352, 3.2258, 2.4248, 1.4265, 1.7188, 3.2992, 4.6512, 3.0395]
    #lambdas_itf = [0.7868, 1.4294, 1.7346, 1.1045, 1.4680, 3.6751, 0.6734, 0.6219, 0.6262, 0.5896]
    # Single null only (no A31 and A32).
    probes      = ['A2',  'A17',  'A18',  'A19',  'A28', 'A33',  'A34',  'A35']
    lambdas_otf = [1.2547, 4.1982, 3.2352, 3.2258, 2.4248, 3.2992, 4.6512, 3.0395]
    lambdas_itf = [0.7868, 1.4294, 1.7346, 1.1045, 1.4680, 0.6219, 0.6262, 0.5896]

    bar_width = 0.35
    fontsize = 24
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    index = np.arange(0, len(probes))
    ax1.bar(index, lambdas_otf, bar_width, color='b', label='OTF', alpha=0.7)
    ax1.bar(index+bar_width, lambdas_itf, bar_width, color='r', label='ITF', alpha=0.7)
    ax1.set_xticks(index + bar_width/2)
    ax1.set_xticklabels(probes, fontsize=fontsize, weight='bold')
    ax1.set_ylabel('Radial Falloff Length (cm)', fontsize=fontsize, weight='bold')
    ax1.legend(fontsize=fontsize, frameon=True, edgecolor='k')
    ax1.tick_params(axis='y', labelsize=fontsize)
    fig.tight_layout()
    fig.show()
