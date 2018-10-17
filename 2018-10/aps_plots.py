import pretty_plots as pp
import numpy as np
import pandas as pd


def itf_otf_plot():
    lambda_ne = np.array([4.605, 4.544, 3.410, 5.623, 3.058, 2.164, 1.722, 2.612,
                          2.401, 2.027, 1.453])
    lambda_ne_err = np.array([0.4552, 0.5542, 0.7938, 1.2043, 0.4252, 0.3956,
                              0.2197, 0.7343, 0.4476, 0.2407, 0.251])
    itf_otf = np.array([3.323809525, 2.302325577, 1.733333322, 2.91666666,
                        0.749999985, 1.850000003, 0.86734694, 0.800000011,
                        0.727272725, 1.375000004, 0.345238096])
    itf_otf_err = np.array([0.393703608, 0.462717516, 0.445963232, 0.6249957,
                            0.406947537, 0.589904423, 0.138455898, 0.372362303,
                            0.222509805, 0.420682602, 0.083429511])

    fig = pp.pplot(x=lambda_ne, xerr=lambda_ne_err, y=itf_otf, yerr=itf_otf_err,
                   xlabel='Density Falloff Length (cm)', ylabel='ITF/OTF Max W Content')


def twod_avg_plots():
    filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/My Slides and Sheets/2018-10/avg_lams_profs.xlsx'
    df = pd.read_excel(filename, sheet_name=3)

    avgc_x = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5])
    avgc_y = np.array([0.823235462, 0.521010541, 0.474225485, 0.49210778,
                       0.466843473, 0.459609305, 0.472382558, 0.500823746,
                       0.581662354, 0.613154856, 0.831188867])

    avgb_x = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5,
                       2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5])
    avgb_y = np.array([0.838255109, 0.680718111, 0.653292943, 0.633479284,
                       0.615370476, 0.600568093, 0.587378014, 0.595564799,
                       0.583085049, 0.572468541, 0.534090827, 0.546292564,
                       0.553421959, 0.558379920, 0.561794678, 0.596267804,
                       0.562282442, 0.564942857, 0.578859443, 0.593684058,
                       0.742640932])

    test68_x = np.array([0.2499999, 0.75, 1.25, 2.5, 3.75, 4.25, 4.7500001])
    test68_y = np.array([1, 0.972940722, 0.895841523, 0.840417983, 0.902008185,
                         0.965122997, 0.993934531])


    figb = None
    figc = None
    for col in df.columns:
        if col[0] == 'b' and col[-1] == 'x':
            pname = col[:3]
            x = df[pname+'_x'].values
            y = df[pname+'_y'].values
            figb = pp.pplot(x, y, fmt='k-', show_fig=False,
                            fig=figb, lw=3, alpha=0.25)

        elif col[0] == 'c' and col[-1] == 'x':
            pname = col[:3]
            x = df[pname+'_x'].values
            y = df[pname+'_y'].values
            figc = pp.pplot(x, y, fmt='k-', show_fig=False,
                            fig=figc, lw=3, alpha=0.25)



    figb = pp.pplot(test68_x, test68_y, label='3DLIM Simulation', fig=figb, fmt='-', color=8)
    figb = pp.pplot(avgb_x, avgb_y, '-', label='Average', fig=figb, xlabel='Z Position (mm)',
                    ylabel='Normalized W Counts')
    figb.axes[0].legend(loc='lower left', fontsize=26)
    figb.tight_layout()
    figb.show()

    figc = pp.pplot(avgc_x, avgc_y, '-', label='Average', fig=figc, xlabel='Z Position (mm)',
                    ylabel='Normalized W Counts')
    figc.axes[0].legend(loc='lower left', fontsize=26)
    figc.tight_layout()
    figc.show()

#itf_otf_plot()
twod_avg_plots()
