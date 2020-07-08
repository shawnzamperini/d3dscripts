import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pretty_plots as pp
import sys
from scipy.interpolate import interp1d


au03_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AU03_Map_Analysis.xlsx'
ad03_file = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/AD03_Map_Analysis.xlsx'

plist = ['A03', 'A04']

for p in plist:

    # Load each U and D LAMS file.
    uname = 'A' + 'U' + p.split('A')[1]
    dname = 'A' + 'D' + p.split('A')[1]
    upath = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/{}_Map_Analysis.xlsx'.format(uname)
    dpath = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
            'Polodial_Scans/New Map Script Results/{}_Map_Analysis.xlsx'.format(dname)
    dfu = pd.read_excel(upath)
    dfd = pd.read_excel(dpath)

    # Do the pivot thing to have the index as axial values.
    dfu = dfu.pivot(columns='z Location [mm]', index='Axial Location [mm]')
    dfd = dfd.pivot(columns='z Location [mm]', index='Axial Location [mm]')

    # Then the LAMS values.
    u_x     = dfu.index.values / 10  # mm to cm
    u_y     = dfu.mean(axis=1).values
    u_y_err = dfu.std(axis=1).values
    d_x     = dfd.index.values / 10  # mm to cm
    d_y     = dfd.mean(axis=1).values
    d_y_err = dfd.std(axis=1).values

    # Assign to correct side.
    if p in ['A03', 'A04', 'A05', 'A06', 'A07']:
        itf_x = d_x
        itf_y = d_y
        itf_y_err = d_y_err
        otf_x = u_x
        otf_y = u_y
        otf_y_err = u_y_err
    elif p in ['A08', 'A33', 'A34', 'A35']:
        itf_x = u_x
        itf_y = u_y
        itf_y_err = u_y_err
        otf_x = d_x
        otf_y = d_y
        otf_y_err = d_y_err

    # Do interpolation functions so we can divide the two sides.
    f_itf     = interp1d(itf_x, itf_y)
    f_itf_err = interp1d(itf_x, itf_y_err)
    f_otf     = interp1d(otf_x, otf_y)
    f_otf_err = interp1d(otf_x, otf_y_err)

    x_int       = np.linspace(0, 5, 100)
    itf_int     = f_itf(x_int)
    itf_err_int = f_itf_err(x_int)
    otf_int     = f_otf(x_int)
    otf_err_int = f_otf_err(x_int)

    # Calculate ITF/OTF.
    itfotf = itf_int / otf_int
    itfotf_err = itfotf * np.sqrt(np.square(itf_err_int/itf_int)
                 + np.square(otf_err_int/otf_int))

    # Print out values.
    print('Distance along probe:')
    for i in x_int:
        print(i)
    print('ITF/OTF')
    for i in itfotf:
        print(i)
    print('ITF/OTF Error')
    for i in itfotf_err:
        print(i)

    # Plot it.
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, figsize=(15,10))
    ax1.plot(x_int, itf_int, 'r', label='ITF')
    ax1.plot(x_int, otf_int, 'b', label='OTF')
    ax2.plot(x_int, itfotf)
    ax1.set_xlabel('Distance along probe (cm)')
    ax1.set_ylabel('LAMS Counts')
    ax2.set_ylabel('ITF/OTF')
    fig.tight_layout()
    fig.show()
