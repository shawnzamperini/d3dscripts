import LimPlots
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Century Gothic"
empty = False

# Just do a simple comparison with a flat source.
ncpath527 = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-020.nc"
ncpath267 = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-006.nc"
ncpath527_flat = "/Users/zamperini/Documents/d3d_work/lim_runs/184527/mcp-184527-021.nc"
ncpath267_flat = "/Users/zamperini/Documents/d3d_work/lim_runs/184267/mcp-184267-007.nc"
if False:
    #absfac = 4.749e14 # For total deposition need to multiply by 2.3, but we aren't doing that here.
    #absfac = 2.330e+14
    absfac = 2.375e+15 * 300
else:
    absfac = 1.0

print("Warning! Swapping 527 out for 267 for PSI plots")
ncpath527 = ncpath267
ncpath527_flat = ncpath267_flat

lp527 = LimPlots.LimPlots(ncpath527)
lp527flat = LimPlots.LimPlots(ncpath527_flat)

# Grab the data. We will do custom plots for our usage.
d527 = lp527.plot_par_rad("nz", 21, charge="all", showplot=False)
d527flat = lp527flat.plot_par_rad("nz", 21, charge="all", showplot=False)

vmin = (np.array((d527["Z"], d527flat["Z"])) * absfac).min()
vmax = (np.array((d527["Z"], d527flat["Z"])) * absfac).max()

# Convert the radial coordinate to R-Rsep.
rorigin = 2.318
rsep = 2.26048
X = d527["X"]
Y = rorigin - d527["Y"] - rsep

negboundx = d527["neg_bound_x"]
posboundx = d527["pos_bound_x"]
boundy = rorigin - d527["bound_y"] - rsep

cmap = "magma"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)

cont1 = ax1.pcolormesh(X, Y, d527flat["Z"]*absfac, shading="auto", vmin=vmin, vmax=vmax, cmap=cmap)
cont2 = ax2.pcolormesh(X, Y, d527["Z"]*absfac, shading="auto", vmin=vmin, vmax=vmax, cmap=cmap)
plt.gca().invert_yaxis()

# We want to fill the areas outside of bounds with grey, overwriting any
# particles that slipped through the cracks.
ax1.fill_betweenx(boundy, -100, negboundx, color="grey", step="pre")
ax1.fill_betweenx(boundy, 100, posboundx, color="grey", step="pre")
ax2.fill_betweenx(boundy, -100, negboundx, color="grey", step="pre")
ax2.fill_betweenx(boundy, 100, posboundx, color="grey", step="pre")

fig.subplots_adjust(right=0.87)
cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
cbar = fig.colorbar(cont1, cax=cbar_ax)
cbar.set_label(r"$\mathdefault{^{13}C\ Density\ (m^{-3})}$", fontsize=14)

ax1.set_xlim([-40, 15])
ax1.set_xlabel("Parallel to B (m)", fontsize=14)
ax2.set_xlabel("Parallel to B (m)", fontsize=14)
ax1.set_ylabel("R-Rsep (m)", fontsize=14)
fig.show()


# Just the not-flat distribution one.
fig, ax2 = plt.subplots(1, 1, figsize=(8, 4))
if not empty:
    cont2 = ax2.pcolormesh(X, Y, d527["Z"]*absfac, shading="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    #cbar = fig.colorbar(cont2, ax=ax2)
    #cbar.set_label(r"$\mathdefault{^{13}C\ Density\ (m^{-3})}$", fontsize=14)
    #cbar.set_label(r"$\mathdefault{^{13}C\ Density\ (Arbitrary)}$", fontsize=14)
plt.gca().invert_yaxis()
ax2.fill_betweenx(boundy, -100, negboundx, color="grey", step="pre")
ax2.fill_betweenx(boundy, 100, posboundx, color="grey", step="pre")
#ax2.set_xlim([-40, 15])
ax2.set_xlim([-17, 40])
ax2.set_ylim([0.1555, 0.029])
ax2.set_xlabel("Parallel to B (m)", fontsize=14)
ax2.set_ylabel("R-Rsep (m)", fontsize=14)
fig.tight_layout()
fig.show()
