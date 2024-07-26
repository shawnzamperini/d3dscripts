import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.colors as colors


# Constants
eps0 = 8.85e-12
ev = 1.602e-19  # C
amu = 1.660e-27  # kg

# Some prefactors for ln(alpha) and nu.
lnalpha_fact = np.log(np.power(eps0, 3/2) * 4 * np.pi * amu / np.power(ev, 5/2))
nu_fact = np.power(ev, 4) * np.power(amu, 3/2) / (4 * np.pi * 
    np.square(eps0) * np.square(amu) * np.power(2 * ev, 3/2))

def calc_nu_lambda_zi(Zz=1, Zi=1, ne=1e19, ni=1e19, Te=10, Ti=10, 
    mz_amu=183.38, mi_amu=2.014):

    # Reduced mass
    mr_amu = mi_amu * mz_amu / (mi_amu + mz_amu)

    # Impurity thermal velocity
    vtz = np.sqrt(2 * Ti * ev / (mz_amu * amu))

    # ln(lambda) factor
    lnalpha = lnalpha_fact + np.log(np.sqrt(Te / ne) * np.square(vtz) * 
        mr_amu / (Zz * Zi))

    # Collision frequency (nu)
    nu = nu_fact * np.square(Zz) * np.square(Zi) * ni * lnalpha / (mz_amu 
        * mr_amu) / (np.power(Ti / mz_amu, 3/2) + 1.3 * np.power(Ti /
        mi_amu, 3/2))

    # W-D colllisional mean free path (lambda).
    lamb = 1 / nu * vtz

    return nu, lamb

# Other constants for these calculations.
ne = 1e19
ni = ne
Te = 10
Ti = Te
mz_amu = 183.38  # W
mi_amu = 2.014   # D
Zi = 1           # D
Zz_max = 74      # W

# Loop through and calculate the W-D collision frequency for each charge state.
Zz_vals = np.arange(1, Zz_max+1)
nus = []
lambs = []
for Zz in Zz_vals:
    nu, lamb = calc_nu_lambda_zi(Zz=Zz, Zi=Zi, ne=ne, ni=ni, Te=Te, Ti=Ti, 
        mz_amu=mz_amu, mi_amu=mi_amu)
    nus.append(nu)
    lambs.append(lamb)

# Some values for the plot
fontsize = 16
lw = 4
conn_length = 10
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(Zz_vals, lambs, lw=lw+1, color="k")

# This is from a matplotlib example for applying a colormap to a line.
points = np.array([Zz_vals, lambs]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

vals_for_norm = np.array(lambs) / conn_length
norm = colors.LogNorm(0.005, vals_for_norm.max())
lc = LineCollection(segments, cmap='RdYlGn_r', norm=norm)

# Set the values used for colormapping
lc.set_array(vals_for_norm)
lc.set_linewidth(lw)
line = ax.add_collection(lc)
#fig.colorbar(line, ax=ax)

# Plot of the results.
ax.axhline(conn_length, linestyle="--", lw=lw, color="k")
#ax.plot(Zz_vals, lambs, lw=lw, color="tab:red")
ax.set_xlabel("W Charge State", fontsize=fontsize)
ax.set_ylabel("W-D Collision MFP (m)", fontsize=fontsize)
ax.set_yscale("log")
ax.set_xlim([0, 20])
ax.set_ylim([3e-3, None])
ax.tick_params(axis='both', which='major', labelsize=fontsize-2)
fig.tight_layout()
fig.show()
