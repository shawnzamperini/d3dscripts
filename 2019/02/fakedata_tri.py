import numpy as np
import pretty_plots as pp


rminrsep = np.linspace(18, 6, 23)
wareal   = np.array([0., 0.001363767, 0.002727533, 0.006818833, 0, 0.0040913,
                     0.0040913, 0.002727533, 0.008182599, 0.008182599, 0.013637666,
                     0.03545793, 0.040912997, 0.07500716, 0.085917293, 0.11455639,
                     0.143195488, 0.175012889, 0.217728965, 0.245, 0.280, 0.325, 0.384])

np.random.seed(19680801)

wareal_low = wareal[:15] + 0.1 * wareal[:15] * np.random.randn(len(wareal[:15]))
wareal_med = wareal[:19] + 0.1 * wareal[:19] * np.random.randn(len(wareal[:19]))
wareal_hig = wareal + 0.1 * wareal * np.random.randn(len(wareal))

#pp.show_colors()
fig = pp.pplot(rminrsep, wareal_hig, fmt='-', label='High Tri.', color=18)
fig = pp.pplot(rminrsep[:19], wareal_med, fmt='-', label='Med Tri.', fig=fig,
               color=8)
fig = pp.pplot(rminrsep[:15], wareal_low, fmt='-', label='Low Tri.', color=6,
               fig=fig, xlabel='R-Rsep (cm)', ylabel='W Areal Density (1e15 cm-2)')
