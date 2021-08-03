import numpy as np
import pretty_plots as pp


# Number of lambda_ne's away
lamb_away_rev = np.array([1.21337012, 1.428339525, 0.979622507, 1.339554147, 1.544795213])
lamb_away_rev_err = lamb_away_rev * 0.1

lamb_away_for = np.array([3.925426957, 4.091525235, 2.096953011, 2.461426518,
                          5.36517793,  2.461133273, 7.822517056, 3.843118169,
                          3.026484426, 1.626114204, 2.300106928, 1.720898512])
lamb_away_for_err = lamb_away_for * 0.1

# The ITF/OTF data
itfotf_rev = np.array([3.77, 5.25, 2.22, 2.50, 2.86])

# A minor fix. The data for A8 (2.86), could reasonably be 2.23 by excluding the
# ITF points that extend past where the OTF stops.
itfotf_rev = np.array([3.77, 5.25, 2.22, 2.50, 2.23])
itfotf_rev_err = itfotf_rev * 0.05  # About right, on the high side for the other reverse probes.

itfotf_for = np.array([0.35, 0.45, 0.41, 0.63, 0.29, 0.87, 0.82, 0.53, 0.66,
                       1.08, 1.15, 1.42])
itfotf_for_err = np.array([0.04, 0.05, 0.05, 0.09, 0.05, 0.07, 0.05, 0.05,
                           0.15, 0.11, 0.12, 0.12])

xlabel = "# of " + r'$\mathrm{\lambda_{ne}}$' + "'s from separatrix"
ms = 10
fig = pp.pplot(lamb_away_rev, itfotf_rev, xerr=lamb_away_rev_err, yerr=itfotf_rev_err, label='Unfavorable', color=8, ms=ms)
fig = pp.pplot(lamb_away_for, itfotf_for, xerr=lamb_away_for_err, yerr=itfotf_for_err, label='Favorable',
               xlabel=xlabel, ylabel='ITF/OTF Total W', fig=fig, logy=False, ms=ms)
