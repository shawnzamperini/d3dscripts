import pretty_plots as pp
import numpy as np


cu4_pol = np.array([2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
cu4_avg = np.array([1, 0.445955458, 0.392902679, 0.373447132, 0.359417752,
                    0.351179181, 0.354183656, 0.342983564, 0.348991165,
                    0.362514404, 0.394384727, 0.796392127])

pp.pplot(cu4_pol, cu4_avg, fmt='-', color=8, xlabel='Poloidal (mm)', ylabel='LAMS Counts (Normalized)')
