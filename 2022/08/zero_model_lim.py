# Same as englehardt_solution just different name with different (hopefully more
# correct) equations.
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Inputs.
r0 = -0.5
rsep = 0.0
rwall = 0.10
rlim = 0.08  # 3 cm wide limiter.
riz = rlim - 0.01  # 1 cm ionization length.
#dperp = 5.0
dperp = [1, 1, 5, 10]
#te = np.array([100, 50, 10, 5])  # Te for region [a, b, c, d]
#te = np.array([10, 10, 10, 10])
#ne = np.array([1e19, 8e18, 5e18, 1e18])
#ne = np.array([1e18, 1e18, 1e18, 1e18])
conns = np.array([100, 50, 10, 5])  # Lconn for region [a, b, c, d]
#conns = np.array([50, 50, 50, 50])  # Lconn for region [a, b, c, d]
neut_flux = 1e20
nr = 50

# Create profiles of ne, Te and Ti based off basic prescriptions.
rprof = np.linspace(rsep, rwall, 100)
nesep = 1e19
lambda_ne = 0.05
neprof = nesep * np.exp(-rprof / lambda_ne)
tesep = 100
lambda_te = 0.05
teprof = tesep * np.exp(-rprof / lambda_te)
timult = 3
tiprof = teprof * timult

# Interpolation functions so we can choose the corresponding ne, Te, Ti values
# for the model.
fne = interp1d(rprof, neprof)
fte = interp1d(rprof, teprof)
fti = interp1d(rprof, tiprof)

# Get values for each region just using the innermost values.
ne = np.array([nesep])
te = np.array([tesep])
ti = np.array([tesep * timult])
for r in [rsep, riz, rlim]:
    ne = np.append(ne, fne(r))
    te = np.append(te, fte(r))
    ti = np.append(ti, fti(r))

# Values from mixed material model, just do placeholders for now.
fz = 0.02

# Constants, one-time calculations.
mi = 931.49e6
cs = np.sqrt((te+ti)/mi) * 3e8

# Root for the equations.
x = np.sqrt(cs / (dperp * conns))

# Region d solution.
rsd = np.linspace(rlim, rwall, nr)
nzd = fz * ne[3] * (np.exp(x[3]*rsd) - np.exp(2*x[3]*rwall) * np.exp(-x[3]*rsd)) \
    / (np.exp(x[3]*rlim) - np.exp(2*x[3]*rwall) * np.exp(-x[3]*rlim))

# Region c solution.
rsc = np.linspace(riz, rlim, nr)
C1 = (fz * ne[3] * np.exp(x[2]*rlim) * np.exp(-x[2]*riz) - neut_flux / dperp[2]) \
    / (x[2]*np.exp(x[2]*riz) + x[2]*np.exp(2*x[2]*rlim) * np.exp(-x[2]*riz))
nzc = C1 * (np.exp(x[2]*rsc) - np.exp(2*x[2]*rlim) * np.exp(-x[2]*rsc)) \
    + fz * ne[3] * np.exp(x[2]*rlim) * np.exp(-x[2]*rsc)

# Region b solution.
rsb = np.linspace(rsep, riz, nr)
C2 = (C1 * (np.exp(x[2]*riz) - np.exp(2*x[2]*rlim) * np.exp(-x[2]*riz)) \
    + fz * ne[3] * np.exp(x[2]*rlim) * np.exp(-x[2]*riz)) \
    / (np.exp(x[1]*riz) + np.exp(2*x[2]*rsep) * np.exp(-x[1]*riz))
nzb = C2 * (np.exp(x[1]*rsb) + np.exp(2*x[1]*rsep) * np.exp(-x[1]*rsb))

# Region a solution.
rsa = np.linspace(r0, rsep)
nza = np.full(nr, nzb[0])

print("Core density: {:.2e} m-3".format(nza[0]))

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
