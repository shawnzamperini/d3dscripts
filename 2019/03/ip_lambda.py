import numpy as np
import pretty_plots as pp


ip = np.array([0.995, 1.472, 1.472, 0.506, 0.508, 0.989, 0.989, 0.987, 1.474,
               0.697, 0.699, 0.989, 0.990, 0.991, 0.994, 0.993, 0.992,
               0.991, 0.988, 0.988])
lambda_ne = np.array([ 6.07,  4.78,  4.95, 15.15, 21.60, 11.26,  8.41,  8.02,
                       6.57,  10.86, 11.34, 8.20,  8.31,  7.81,  9.47,  9.12,
                       11.13, 10.11, 6.44, 7.03]) / 10.0

fig = pp.pplot(ip, lambda_ne, xlabel='IP (MA)', ylabel=r'$\mathrm{\lambda_{ne}}$' + ' (cm)')
