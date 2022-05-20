# This script is a super basic demonstration of how to "solve" the equation:
#       dx/dt = sin(t) + t^3            x(0) = -1
# via a Monte-Carlo typie method. The solution is
#      x(t) = -cos(t) + t^4 / 4
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from sklearn.metrics import mean_squared_error
from tqdm import tqdm
from scipy.interpolate import interp1d


N = 500  # Number of points in each solution.
domain = [0, 10]
NN = 20  # Number of solutions to generate.
smooth = True
width_mult = 5
NNN = 10  # Number of iterations, where each iteration uses the best previous solution.

# Idea is to define a number of bins for our solution, and then sample a number
# of values from an initial guess solution to minimize the RMSE between the
# LHS and RHS of the equation.
ts = np.linspace(domain[0], domain[1], 100)

# Using the ODE, fill out values of the dependent variable.
dxdts = np.sin(ts) + np.power(ts, 3)

for i in range(0, N):
    ts_bins = [[] for t in ts]

    for j in range(0, len(ts)):
        t = ts[j]

        # Fill in the ts_bins with guesses of what x could be.
        if i == 0:
            #guess_x = np.power(t, 4) / 4
            guess_x = np.power(t, 3) / 3
        else:
            guess_x = f_best(t)

        center = guess_x
        width = np.abs(center * width_mult)
        x_ran = np.random.normal(center, width)
        ts_bins[j].append(x_ran)

    # Now, at each t, use the two bounding points to get the point-slope equation
    # of dxdt = mt + b. x at this t location is then x = 1/2 * mt^2 + bt + c
    # with the boundary condition x(0) = -1 to determine c.
    for j in range(0, len(ts)):

        # The first and last points we just have to use two points.
        #if j == 0:
        #    ts_sub = ts[0:2]
        #    dxdts_sub = dxdts[0:2]
        #elif j == N-1:
        #    ts_sub = ts[-2:]
        #    dxdts_sub = dxdts[-2:]
        #else:
        #    ts_sub = ts[j-1:j+2]
        #    dxdts_sub = dxdts[j-1:j+2]
        #t0 = ts[j]

        # Do linear fit, but then what's the BC for this "temporary" integral of
        # the linear fit to determine c? Let's try with an intial guess on what
        # x could be.
        #m, b = np.polyfit(ts_sub, dxdts_sub, 1)

        # For each educated guess "random" x value in this tbin location, identify
        # the one that minimizes the error between the LHS and RHS
        for x_ran in ts_bins[j]:
            #c = x_ran - 0.5 * m * t**2 - b * t
            #x = 0.5 * m * t**2 + b * t + c
        results["x"][i] = x
