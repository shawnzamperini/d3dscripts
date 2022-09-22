# Plot the full connection lengths for a cross section.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


#path = "/Users/zamperini/Documents/d3d_work/mafot_files/170848/lam_full_vessel240.dat"
path = "/Users/zamperini/Documents/d3d_work/mafot_files/190175/lam_lobes.dat"
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
df = pd.read_csv(path, skiprows=52, names=columns, delimiter="\t")

# From input I know that there are 275 R points and 544 Z points equally spaced.
r = np.linspace(0.98, 1.97, 198)
z = np.linspace(0.7, 1.37, 134)
#r = df["R (m)"].values
#z = df["Z (m)"].values
R, Z = np.meshgrid(r, z)
L = df["Lconn (km)"].values.reshape(R.shape)

norm = LogNorm(vmin=0.005, vmax=1)

fig, ax1 = plt.subplots()

ax1.set_title("190175.3500")
ax1.pcolormesh(R, Z, L, norm=norm, cmap="inferno")
ax1.set_aspect("equal")

fig.tight_layout()
fig.show()
