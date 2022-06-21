# Plot the full connection lengths for a cross section.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle


path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/lam_dimes_err.dat"
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
df = pd.read_csv(path, skiprows=52, names=columns, delimiter="\t")

# From input I know that there are 275 R points and 544 Z points equally spaced.
r = np.linspace(1.01, 1.9, 222)
z = np.linspace(-1.365, -0.8, 141)
#r = df["R (m)"].values
#z = df["Z (m)"].values
R, Z = np.meshgrid(r, z)
L = df["Lconn (km)"].values.reshape(R.shape)

norm = LogNorm(vmin=0.005, vmax=0.5)

# Rectangles for metal ring and DiMES.
dimes = Rectangle((1.466, -1.25-0.02), 1.514-1.466, 0.02, color="k")

fig, ax1 = plt.subplots()

ax1.pcolormesh(R, Z, L, norm=norm, cmap="inferno")
ax1.set_aspect("equal")
ax1.add_patch(dimes)

fig.tight_layout()
fig.show()
