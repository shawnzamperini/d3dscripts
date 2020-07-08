import pretty_plots as pp
import numpy as np


peaking_f = np.array([0.413,       0.734,       0.790,      0.929,       1.036,       0.564,       0.768,       0.984])
total_f   = np.array([0.303030303, 0.632727273, 0.52345679, 0.689570552, 0.44891945,  0.355293105, 0.515093095, 0.85257145])
peaking_r = np.array([1.072,      1.478,       1.670,       0.994])
total_r   = np.array([1.06741573, 1.290322581, 1.318587466, 1.592779647])
peaking_3d = np.array([1.36, 1.07, 0.99, 0.94, 0.72])
total_3d   = np.array([1.91, 1.08, 1.01, 0.94, 0.5])
error = 0.05

fig = pp.pplot(peaking_f, total_f, xerr = error * peaking_f, yerr = total_f * error, color=6, label='Favorable', ms=12)
fig = pp.pplot(peaking_r, total_r, xerr = error * peaking_r, yerr = total_r * error, color=8, fig=fig, label='Unfavorable', ms=12)
fig = pp.pplot(peaking_3d, total_3d, fmt='*', color=18, ms=22, fig=fig, xlabel='OTF/ITF Peaking', ylabel='ITF/OTF Total', label='3DLIM')

fig.axes[0].axhline(1.0, linestyle='--', color='k', lw=2, alpha=0.75)
fig.axes[0].axvline(1.0, linestyle='--', color='k', lw=2, alpha=0.75)
