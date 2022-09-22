# This script tries to solve the ODE's in Dennis's paper. It uses the following
# cookbook example as a guide:
# https://scipy-cookbook.readthedocs.io/items/CoupledSpringMassSystem.html
#
# Useful StackOverflow:
# https://stackoverflow.com/questions/65344347/coupled-system-of-4-differential-equations-python

from scipy.integrate import odeint
import numpy as np
import sympy as sm
from sympy import symbols
import matplotlib.pyplot as plt


# The units are in in normalized units. E.g.,
# s = s / sO
# n = n / nO
# A = A / AO
# S = S / SO
# M = doesn't need it
# Where O indicates the measurement location. Dennis assumes this is the
# stagnation point, but we do not have to assume that.

# The first issue is that we need to decouple the equations. As written, dn/ds
# is dependent on dM/ds, which means as it stands we can't use. We use sympy
# to solve for dM/ds. This only needs to be run once, but I include it here
# for completeness.

# Define symbols.
#dMds, M, S, n, cs, A, dAds, dnds = symbols("dMds, M, S, n, cs, A, dAds, dnds", real=True)

# Declare RHS as expressions.
#M1 = (1-M**2)**-1 * (S * (1+M**2) / (n*cs) - M * (A**-1 * dAds))
#n1 = (1+M**2)**-1 * (-2*n*M*dMds - M**2*n*(A**-1 * dAds))

# Solver looks for equations of the form, e.g., dM/ds - M1 = 0.
#result = sm.solve([dMds - M1, dnds - n1], [dMds, dnds])
#print(result)

# Define constants.
Te = 7  # eV
#Te = 1
rt_rO = 0.75  # At / AO = Rt / RO, this define the linear shape of A.
m = 931.49e6 / 3e8**2  # eV s2 m-2
#m = 1
num_points = 100
#S = 1
S = 0

# ODE solver parameters.
abserr = 1.0e-8
relerr = 1.0e-6

# Setup up needed variables.
cs = np.sqrt(2 * Te / m)   # Sound speed.
#cs = 1
s = np.linspace(0, 1.0, num_points)  # The domain.
#dAds = (rt_rO - 1) / (1 - 0)   # Slope for A.
dAds = 0
A = dAds * s - dAds * 0 + 1  # Values for A across the domain.
#A = 1.0

def getA(sval):
    idx = np.argmin(np.abs(s-sval))
    #return A[idx]
    return 1.0

def vectorfield(w, s, p):
    """
    Defines the differential equations for the coupled M, n equations.

    Arguments:
        w : vector of the state variables:
            w = [M, n]
        s : Normalized distance
        p : Vector of the parameters:
            p = [A, cs, S]
    """

    M, n = w
    A, cs, S = p

    n2 = 1/n

    # Create f = (M', n'). Equations are actually defined by the above sympy solver.
    #f = [(-A*M**2*S - A*S + M*cs*dAds*n)/(A*M**2*cs*n - A*cs*n),
    #     (2*A*M*S - M**2*cs*dAds*n)/(A*M**2*cs - A*cs)]

    f = [(-getA(s)*M**2*S*cs**2*m - 2*getA(s)*S*Te + 2*M*Te*cs*dAds*n)/(getA(s)*M**2*cs**3*m*n - 2*getA(s)*Te*cs*n),  # dM/ds
         #(-A*M**2*S*cs**2*m*n2 - 2*A*S*Te*n2 + 2*M*Te*cs*dAds/n2)/(A*M**2*cs**3*m - 2*A*Te*cs),  # n = 1/n2
         (2*getA(s)*M*S*cs*m - M**2*cs**2*dAds*m*n)/(getA(s)*M**2*cs**2*m - 2*getA(s)*Te)]  # dn/ds
         #(2*A*M*S*cs*m - M**2*cs**2*dAds*m/n2)/(A*M**2*cs**2*m - 2*A*Te)]  # n = 1/n2

    #f = [(1+M**2)/(n-M**2), 2*M/(M**2-1)]
    return f

# Initial conditions.
M0 = 0  # Assuming O == stagnation point
n0 = 1e18

# Pack up the parameters and call the solver.
w0 = [M0, n0]
p = [A, cs, S]
wsol = odeint(vectorfield, w0, s, args=(p,))

fig, ax1 = plt.subplots(figsize=(5, 4))

#ax11 = ax1.twinx()
ax1.plot(s, wsol[:,0], label="Mach", color="tab:red", lw=3)
ax1.plot(s, wsol[:,1], label="n/nO", color="tab:green", lw=3)
ax1.legend()
ax1.set_ylabel("Normalized Unit", fontsize=16)
ax1.set_xlabel("Normalized S Coordinate", fontsize=16)
ax1.grid()
ax1.set_xlim([0, 1])
ax1.set_ylim([0, 1])

fig.tight_layout()
fig.show()
