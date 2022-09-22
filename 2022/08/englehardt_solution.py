# This script is a starter sort of rough draft for the analytical solution of
# an Englehardt type model to describe impurity densities when the only source
# is from the wall/limiters.
import numpy as np
import matplotlib.pyplot as plt


# Consider DIII-D-like parameters.
r0 = -0.5
rsep = 0.0
rwall = 0.10
rlim = 0.08  # 3 cm wide limiter.
riz = rlim - 0.01  # 1 cm ionization length.
dperp = 5.0
#te = np.array([100, 10, 5, 5])  # Te for region [a, b, c, d]
te = np.array([10, 10, 10, 10])
#ne = np.array([1e19, 8e18, 5e18, 1e18])
ne = np.array([1e18, 1e18, 1e18, 1e18])
#conns = np.array([100, 50, 10, 5])  # Lconn for region [a, b, c, d]
conns = np.array([50, 50, 50, 50])  # Lconn for region [a, b, c, d]
neut_flux = 1e18

# Values from mixedmaterial model, just do placeholders for now.
fz = 0.02

# Constants, one-time calculations.
mi = 931.49e6
cs = np.sqrt(2 * te / mi) * 3e8
Spn = neut_flux / (rlim - riz)

# The roots in the equations.
x1 = np.sqrt(dperp * cs / conns) / dperp  # Positive root.
x2 = -x1  # Negative root.


# Region (d) solution.
nr = 50
C1 = fz * ne[3] / (np.exp(x1[3]*rlim) - np.exp((x1[3]-x2[3])*rwall) * np.exp(x2[3]*rlim))
rsd = np.linspace(rlim, rwall, nr)
nzd = C1 * (np.exp(x1[3]*rsd) - np.exp((x1[3]-x2[3])*rwall) * np.exp(x2[3]*rsd))

# Region (c) solution.
#B = ((fz * ne[3] * np.exp(-x2[2]*rlim) - Spn*conns[2]/cs[2] * np.exp(-x2[2]*rlim)) \
#    * np.exp(x2[2]*riz * np.exp(-x2[1]*riz)) + (np.exp(-x2[1]*riz) - 1) * Spn*conns[2]/cs[2]) \
#    * (np.exp(x1[1]*riz) - np.exp((x1[1]-x2[1])*rlim) * np.exp(x2[1]*riz) \
#    - np.exp(x1[2]*riz) * np.exp(-x2[1]*riz) + np.exp((x1[2]-x2[2])*rlim) \
#    * np.exp(x2[2]*riz) * np.exp(-x2[1]*riz))**(-1) * (np.exp(x1[2]*riz) \
#    - np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*riz)) + Spn*conns[2]/cs[2] \
#    + (fz * ne[3] * np.exp(-x2[2]*rlim) - Spn*conns[2]/cs[2] * np.exp(-x2[2]*rlim)) \
#    * np.exp(x2[2]*riz)

#C2 = (B - Spn * conns[2] / cs[2] - (fz * ne[3] * np.exp(-x2[2]*rlim) \
#    - Spn * conns[2] / cs[2] * np.exp(-x2[2]*rlim)) * np.exp(x2[2]*riz)) / \
#    (np.exp(x1[2]*riz) - np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*riz))

C2 = (-neut_flux / dperp - x2[2]*fz * ne[3] * np.exp(-x2[2]*rlim) * np.exp(x2[2]*riz)) \
    / (x1[2] * np.exp(x1[2]*riz) - x2[2] * np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*riz))

rsc = np.linspace(riz, rlim, nr)
#nzc = C2 * (np.exp(x1[2]*rsc) - np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*rsc)) + \
#    (fz * ne[3] * np.exp(-x2[2]*rlim) - Spn * conns[2] / cs[2] \
#    * np.exp(-x2[2]*rlim)) * np.exp(x2[2]*rsc) + Spn * conns[2] / cs[2]
nzc = C2 * (np.exp(x1[2]*rsc) - np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*rsc)) + \
    fz * ne[3] * np.exp(-x2[2]*rlim) * np.exp(x2[2]*rsc)

# Region (b) solution.
#C3 = C2 * (np.exp(x1[2]*riz) - np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*riz)) \
#    * np.exp(-x2[1]*riz) + (fz*ne[3] - Spn*conns[2]/cs[2]) * np.exp(-x2[2]*rlim) \
#    * np.exp(x2[2]*riz) * np.exp(-x2[1]*riz) + Spn*conns[2]/cs[2] * np.exp(-x2[1]*riz)

#C5 = x2[1] * C3 * np.exp(x2[1]*rsep) / (x1[1]*np.exp(x1[1]*rsep) \
#    - x2[1]*np.exp((x1[1]-x2[1])*riz) * np.exp(x2[1]*rsep))

#C4 = (A - C3 * np.exp(x2[1]*rsep)) / (np.exp(x1[1]*rsep) - np.exp((x1[1]-x2[1])*riz) \
#    * np.exp(x2[1]*rsep))

A = (C2 * (np.exp(x1[2]*riz) - np.exp((x1[2]-x2[2])*rlim) * np.exp(x2[2]*riz)) \
    + fz*ne[3]*np.exp(-x2[2]*rlim) * np.exp(x2[2]*riz) + neut_flux / (dperp * x2[1])) \
    * (np.exp(x1[1]*riz) - x1[1] * np.exp(x1[1]*riz) / x2[1])**(-1) \
    * (np.exp(x1[1]*rsep) - x1[1] * np.exp(x1[1]*riz) / (x2[1] * np.exp(x2[1]*riz)) \
    * np.exp(x2[1]*rsep)) - neut_flux / (dperp * x2[1] * np.exp(x2[1]*riz)) * np.exp(x2[1]*rsep)

C4 = (A + neut_flux / (dperp * x2[1] * np.exp(x2[1]*riz))) / (np.exp(x1[1]*rsep) \
    - x1[1] * np.exp(x1[1]*riz) / (x2[1] * np.exp(x2[1]*riz)) * np.exp(x2[1]*rsep))

rsb = np.linspace(rsep, riz, nr)
#nzb = C5 * (np.exp(x1[1]*rsb) - np.exp((x1[1]-x2[1])*riz) * np.exp(x2[1]*rsb)) \
#    + C3 * np.exp(x2[1] * rsb)
nzb = C4 * (np.exp(x1[1]*rsb) - x1[1]*np.exp(x1[1]*riz) / (x2[1]*np.exp(x2[1]*riz)) \
    * np.exp(x2[1]*rsb)) - neut_flux / (dperp * x2[1] * np.exp(x2[1]*riz)) * np.exp(x2[1]*rsb)

# Region (a)
#A = C5 * (np.exp(x1[1]*rsep) - np.exp((x1[1]-x2[1])*riz) * np.exp(x2[1]*rsep)) \
#    + C3 * np.exp(x2[1]*rsep)
rsa = np.linspace(r0, rsep, 2)
nza = np.full(len(rsa), A)



fig, ax = plt.subplots()

ax.axvline(rsep, color="k", linestyle="--")
ax.axvline(riz, color="k", linestyle="--")
ax.axvline(rlim, color="k", linestyle="--")
ax.axvline(rwall, color="k", linestyle="--")
ax.plot(rsd, nzd)
ax.plot(rsc, nzc)
ax.plot(rsb, nzb)
ax.plot(rsa, nza)
ax.set_xlim(-0.02, rwall)

ax.set_xlabel("R-Rsep (m)")
ax.set_ylabel("nz (m-3)")

fig.tight_layout()
fig.show()
