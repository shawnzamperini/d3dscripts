import LimPlots
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LogNorm, SymLogNorm, BoundaryNorm
from matplotlib.lines import Line2D


# Load, pull out data for plots.
# No suffix - Te = 5 Ti = 15
# 002 - Te = 10 Ti = 30
ncpath_sic_c = "/Users/zamperini/Documents/d3d_work/lim_runs/mm-sic-c.nc"
ncpath_c = "/Users/zamperini/Documents/d3d_work/lim_runs/mm-c-only.nc"
r_origin = 0.12   # About 12 cm from the separatrix in 187111
concentration = True

lp_sic = LimPlots.LimPlots(ncpath_sic_c)
lp_c = LimPlots.LimPlots(ncpath_c)
sic = lp_sic.plot_par_rad("nz", 21, showplot=False)
c = lp_c.plot_par_rad("nz", 21, showplot=False)
X = c["X"]
Y = c["Y"]
Z_c = c["Z"]
Z_sic = sic["Z"]

# Pull out ne as well for concentration.
if concentration:
    ne = lp_c.plot_par_rad("ne", 21, showplot=False)

# Need some additional stuff from the NetCDF file.
nc = lp_c.nc
qxs = nc.variables["QXS"][:].data
qedges1 = nc.variables["QEDGES"][:].data[0]
qedges2 = nc.variables["QEDGES"][:].data[1]
nonzero = np.nonzero(qedges1)
qxs = qxs[:len(qedges1)][nonzero]
qedges1 = qedges1[nonzero]
qedges2 = qedges2[nonzero]
qxs = np.append(qxs, 0)
qedges1 = np.append(qedges1, 0)
qedges2 = np.append(qedges2, 0)
cl = float(nc['CL'][:].data)
ca  = float(nc['CA'][:].data)
caw = float(nc['CAW'][:].data)

# Convert coordinates to a rough R-Rsep, and then convert to cm.
Y = r_origin - Y
qxs = r_origin - qxs
X = X * 100
Y = Y * 100
qxs = qxs * 100
qedges1 = qedges1 * 100
qedges2 = qedges2 * 100

# Setup for plotting.

xlim1 = [-0.2*100, 0.2*100]; equal = True
ncolors = 10
cmap = plt.cm.get_cmap('inferno', ncolors)

if concentration:
    norm = LogNorm(vmin=1e-3, vmax=1e0)
    bounds = np.geomspace(1e-3, 1e0, ncolors)
    yticks = np.power(10.0, np.arange(-3, 1e0))
    yticklabels = [f'$10^{{{np.log10(b):.0f}}}$' for b in yticks]
    Z_c = Z_c / ne["Z"]
    Z_sic = Z_sic / ne["Z"]
else:
    norm = LogNorm(vmin=1e15, vmax=1e18)
    #bounds = np.power(10.0, np.arange(14, 19))
    bounds = np.geomspace(1e15, 1e18, ncolors)
    yticks = np.power(10.0, np.arange(15, 19))
    yticklabels = [f'$10^{{{np.log10(b):.0f}}}$' for b in yticks]

norm = BoundaryNorm(boundaries=bounds, ncolors=ncolors)

ax_yticks = np.arange(0.10, 0.18, 0.02) * 100

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.set_xlabel("Distance along field line (cm)", fontsize=14)
ax.set_ylabel("R-Rsep (cm)\n", fontsize=14)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,3))

cmesh = ax1.pcolormesh(X, Y, Z_c, cmap=cmap, norm=norm)
ax1.fill_between(-qedges1, qxs, Y.max(), color="grey")
ax1.fill_between(qedges2, qxs, Y.max(), color="grey")
ax1.plot(-qedges1, qxs, color="grey", lw=2)
ax1.plot(qedges2, qxs, color="grey", lw=2)
ax1.set_xlim(xlim1)
ax1.set_ylim([Y.max(), (r_origin - ca)*100])
ax1.set_aspect("equal")
ax1.set_yticks(ax_yticks)

cmesh = ax2.pcolormesh(X, Y, Z_sic, cmap=cmap, norm=norm)
ax2.fill_between(-qedges1, qxs, Y.max(), color="grey")
ax2.fill_between(qedges2, qxs, Y.max(), color="grey")
ax2.plot(-qedges1, qxs, color="grey", lw=2)
ax2.plot(qedges2, qxs, color="grey", lw=2)
ax2.set_xlim(xlim1)
ax2.set_ylim([Y.max(), (r_origin - ca)*100])
ax2.set_aspect("equal")
ax2.set_yticks(ax_yticks)

# https://stackoverflow.com/questions/13784201/how-to-have-one-colorbar-for-all-subplots
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(cmesh, cax=cbar_ax)
if concentration:
    cbar_ax.set_ylabel(r"$\mathdefault{n_C}$ / $\mathdefault{n_e}$", fontsize=14)
else:
    cbar_ax.set_ylabel("C Density (m-3)", fontsize=14)
cbar_ax.set_yticks(yticks)
cbar_ax.set_yticklabels(yticklabels)

ax1.text(0.85, 0.75, "Graphite", transform=ax1.transAxes, fontsize=14, bbox={"fc":"white"}, ha="center")
ax2.text(0.85, 0.75, "SiC", transform=ax2.transAxes, fontsize=14, bbox={"fc":"white"}, ha="center")

ax1.set_xticks([])

fig.show()
