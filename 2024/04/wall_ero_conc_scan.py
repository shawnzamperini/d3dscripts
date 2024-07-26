import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# Pentration factor.
fp = 0.001

# Line-averaged density.
ne = 1e20

# Machine minor radius.
r = 2.8

# Erosion flux at the wall.
gamma_ero_wall = np.geomspace(1e15, 1e25, 49)

# Core residence time. 
t_res = np.geomspace(0.01, 10, 50)

# Create 2D array of concentration values and calculate.
cz = np.zeros((len(gamma_ero_wall), len(t_res)))
for i in range(cz.shape[0]):
    for j in range(cz.shape[1]):
        cz[i, j] = 2 * fp * gamma_ero_wall[i] * t_res[j] / (r * ne)

# meshgrid for plotting 2D plot. 
X, Y = np.meshgrid(gamma_ero_wall, t_res)

# Critical value to plot a line for.
cz_crit = 1e-5

fontsize = 14
fig, ax = plt.subplots(figsize=(6, 4))
mesh = ax.pcolormesh(X, Y, cz.T, shading="nearest", 
    norm=colors.LogNorm(vmin=1e-8, vmax=1.0), cmap="inferno")
ax.contour(X, Y, cz.T, levels=[cz_crit], colors="k", linewidths=4)
cont = ax.contour(X, Y, cz.T, levels=[cz_crit], colors="r", linewidths=3)
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label("Core Impurity Concentration", fontsize=fontsize)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"Wall Erosion Flux ($\mathdefault{m^{-2} s^{-1}}$)", 
    fontsize=fontsize)
ax.set_ylabel("Impurity Core Residence Time (s)", fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=fontsize-4)
cbar.ax.tick_params(axis='both', which='major', labelsize=fontsize-4)
fig.tight_layout()
fig.show()
