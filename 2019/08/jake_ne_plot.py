import pretty_plots as pp
import numpy as np


nebar_unf = np.array([3.20E+13, 3.12E+13, 3.28E+13, 3.00E+13, 5.19E+13])
nebar_fav = np.array([2.27E+13, 2.37E+13, 5.75E+13, 4.84E+13, 5.01E+13,
                      5.58E+13, 5.73E+13, 7.09E+13, 4.38E+13, 4.88E+13,
                      4.79E+13, 4.63E+13])

neerr_unf = np.array([0.06, 0.10, 0.10, 0.09, 0.10]) * 1E13
neerr_fav = np.array([0.04, 0.08, 0.23, 0.41, 0.15,
                      0.15, 0.15, 0.24, 0.79, 0.13,
                      0.16, 0.16]) * 1E13

#neerr_unf = nebar_unf * 0.1
#neerr_fav = nebar_fav * 0.1

itfotf_unf = np.array([3.77, 5.25, 2.22, 2.50, 2.86])
itfotf_fav = np.array([0.35, 0.45, 0.41, 0.63, 0.29,
                       0.87, 0.82, 0.53, 0.66, 1.08,
                       1.15, 1.42])

yerrs_unf = itfotf_unf * 0.1
yerrs_fav = itfotf_fav * 0.1

fig = pp.pplot(nebar_unf, itfotf_unf, yerr=yerrs_unf, color=8, logy=True, fmt='s', ms=13, label='Bx' + r'$\nabla$' + 'B' +r'$\uparrow$')
fig = pp.pplot(nebar_fav, itfotf_fav, yerr=yerrs_fav, color=6, logy=True, ms=13, label='Bx' + r'$\nabla$' + 'B' +r'$\downarrow$',
               fig=fig, xlabel='Plasma Density (cm-3)', ylabel='ITF/OTF W Content', yrange=[0.1, 10],
               xrange=[1e13, 9e13])

# Just different labels.
fig = pp.pplot(nebar_unf, itfotf_unf, xerr=neerr_unf, yerr=yerrs_unf, color=8, logy=True, fmt='s', ms=13, label='Bx' + r'$\nabla$' + 'B' +r'$\uparrow$')
fig = pp.pplot(nebar_fav, itfotf_fav, xerr=neerr_fav, yerr=yerrs_fav, color=6, logy=True, ms=13, label='Bx' + r'$\nabla$' + 'B' +r'$\downarrow$',
               fig=fig, xlabel='Plasma Density (cm-3)', ylabel='ITF/OTF Total W', yrange=[0.1, 10],
               xrange=[1e13, 9e13])

fig.axes[0].axhline(1.0, color='k', linestyle='--', linewidth=5, alpha=0.7)
