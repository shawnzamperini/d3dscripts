# Script to plot connections lengths for a pari of MAFOT runs looking at the
# window made by the toroidal limiters.
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Rectangle


# Values are at an (R, Z) range of R = 2.365 and Z = [-0.35, .35] for every toroidal angle.
#tor_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_lim_window_lim.dat"
tor_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_lim_window_wide.dat"
pol_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_lim_window_pol.dat"
#pol_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_lim_window_lim.dat"

columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
tor = pd.read_csv(tor_path, skiprows=52, names=columns, delimiter="\t")
pol = pd.read_csv(pol_path, skiprows=52, names=columns, delimiter="\t")

deg = np.arange(0, 360)
deg = deg[::-1]
def get_dat(df, nz=50):
    Z = df["Z (m)"].values.reshape(360, nz)
    L = df["Lconn (km)"].values.reshape(360, nz) * 1000
    L = np.ma.masked_where(L<=0, L)
    return Z, L

tor_Z, tor_L = get_dat(tor)
pol_Z, pol_L = get_dat(pol)
deg_plot = np.full((50, 360), deg).T

vmin = 0.0
vmax = max(tor_L.max(), pol_L.max())

ncolors = 15
bounds = np.linspace(0, vmax, ncolors)
cmap = plt.cm.get_cmap('inferno', ncolors)
norm = BoundaryNorm(boundaries=bounds, ncolors=ncolors)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

cont1 = ax1.contourf(deg_plot, pol_Z, pol_L, cmap=cmap, norm=norm)
cont2 = ax2.contourf(deg_plot, tor_Z, tor_L, cmap=cmap, norm=norm)

ax1.axvline(180, color="skyblue", linestyle="--", lw=3)
ax2.axvline(180, color="skyblue", linestyle="--", lw=3)

ax1.text(150, -0.3, "WITS", rotation="vertical", color="skyblue", fontsize=20)
ax2.text(150, -0.3, "WITS", rotation="vertical", color="skyblue", fontsize=20)

ax1.set_facecolor("grey")
ax2.set_facecolor("grey")

cbar = fig.colorbar(cont1, ax=ax2)
cbar.set_label("Connection Length (m)")

ax1.set_title("Poloidal Bumpers")
ax2.set_title("Toroidal Limiters")

ax1.set_xlabel("Degree")
ax2.set_xlabel("Degree")
ax1.set_ylabel("Z (m)")

fig.tight_layout()
fig.show()
