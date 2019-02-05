import sys
sys.path.append('/home/shawn/d3dscripts/2018-3')

import numpy             as np
import matplotlib.pyplot as plt
import C_on_CW           as sp
from scipy.interpolate import interpolate
from scipy.signal      import savgol_filter

# Load in all the dataframes from the C_on_CW model.
y_df      = sp.yields()
fluxes_df = sp.get_fluxes(tungsten_frac=np.array([0,0,0,0,0,0]))
sput_df   = sp.get_sput_flux(y_df, fluxes_df)
net_df    = sp.calc_net_fluence()

def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c

def myplot(probe='AD'):
    x_net  = net_df[probe + ' R-Rsep (cm)'][6:-1]
    y_net  = net_df[probe + ' Net Fluence (m-2)'][6:-1]
    x_sput = sput_df['R-Rsep (cm)']
    y_sput = sput_df['Sputt. Fluence from C on 90C10W']

    # Smooth the sputtered data before fitting it.
    y_sput_smooth = savgol_filter(y_sput, 11, 3)

    # Want to have sputtered values at the x values of the net data.
    f_sput = interpolate.interp1d(x_sput, y_sput_smooth)
    y_sput_common = f_sput(x_net)

    # Then add net and sputt to get total.
    y_tot = y_net + y_sput_common

    # For the exponential fit, use corresponding parameters.
    if probe == 'AD':
        a = 2.13e7; b = 1.08; c = 17.2
    elif probe == 'AU':
        a = 1.35e8; b = 1.26; c = -6.93
    else:
        print("Bad probe entry.")

    y_fit = exp_fit(x_net, a, b, c) * 10**15

    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.semilogy(x_net, y_net, 'gx', label='Net')
    ax1.semilogy(x_net, y_tot, 'bx', label='Total')
    ax1.semilogy(x_net, y_sput_common, 'rx', label='Sputtered')
    ax1.semilogy(x_net, y_fit, 'g-', label='Exp. Fit')
    ax1.semilogy(x_sput, y_sput_smooth, 'r.')
    ax1.legend()
    ax1.set_xlabel('R-Rsep (cm)')
    ax1.set_ylabel('W Fluence (m-2)')
    ax1.set_title(probe + ' Fluences')
