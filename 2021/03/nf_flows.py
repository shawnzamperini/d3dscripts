# This plot shows parallel Mach numbers to demonstrate qualitatively similar
# flow patterns as that observed in tokamaks.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Get case representative of each direction.
unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
fav_noflow_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-001.nc"
fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-002c.nc"
unf = oedge_plots.OedgePlots(unf_path)
fav_noflow = oedge_plots.OedgePlots(fav_noflow_path)
fav = oedge_plots.OedgePlots(fav_path)

# Pull out along ring profiles of the Mach number.
unf_ring = 59; fav_ring = 36
unf_s, unf_mach = unf.along_ring(unf_ring, "Mach", plot_it=False)
fav_s, fav_noflow_mach = fav_noflow.along_ring(fav_ring, "Mach", plot_it=False)
fav_s, fav_mach = fav.along_ring(fav_ring, "Mach", plot_it=False)
unf_mask = unf_s != 0.0
fav_mask = fav_s != 0.0
unf_s = unf_s[unf_mask]
fav_s = fav_s[fav_mask]
unf_mach = unf_mach[unf_mask]
fav_noflow_mach = fav_noflow_mach[fav_mask]
fav_mach = fav_mach[fav_mask]

# Distance from separatrix at midplane.
unf_mid_dist = unf.nc.variables["MIDIST"][1][unf_ring]
unf_mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(unf_mid_dist*100)
fav_mid_dist = fav.nc.variables["MIDIST"][1][fav_ring]
fav_mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(fav_mid_dist*100)
unf_psin = unf.nc.variables["PSIFL"][:][unf_ring][0]
unf_psin_str = r"$\mathrm{\psi_N}$" + " = {:.2f}".format(unf_psin)
fav_psin = unf.nc.variables["PSIFL"][:][fav_ring][0]
fav_psin_str = r"$\mathrm{\psi_N}$" + " = {:.2f}".format(fav_psin)
print("Unfavorable: {}  psin = {}".format(unf_mid_str, unf_psin))
print("Favorable:   {}  psin = {}".format(fav_mid_str, fav_psin))

# Plot each.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5

fig, ax = plt.subplots(figsize=(5,4))
ax.axhline(0.0, color="k", linestyle="-", lw=lw)
#ax.plot(unf_s, unf_mach, label="Unfavorable", lw=lw, color=colors[2])
ax.plot(fav_s, fav_noflow_mach, label="Favorable", lw=lw, color=colors[3])
ax.plot(fav_s, fav_mach, lw=lw, color=colors[3], linestyle="--")
ax.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax.set_ylabel("Mach number", fontsize=fontsize)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid()
#ax.legend(fontsize=fontsize)
ax.set_ylim([-1.0, 0.5])
#ax.text(0.45, 0.15, unf_mid_str, bbox=dict(color="white"), transform=ax.transAxes, fontsize=fontsize)
#ax.text(0.45, 0.25, unf_psin_str, bbox=dict(color="white"), transform=ax.transAxes, fontsize=fontsize)
#ax.arrow(22.9, 0.02, 0, -(0.054-(-0.15)), head_width=0.05, head_length=0.1, fc="k", ec="k")
#ax.arrow(0.3, 0.53, 0.0, -0.12, transform=ax.transAxes, width=0.01, color="k", zorder=50)
#ax.text(0.35, 0.38, "M = 0.3", bbox=dict(color="white"), transform=ax.transAxes, fontsize=fontsize)
ax.text(0.05, 0.8, "Towards Outer Target", bbox=dict(color="white"), transform=ax.transAxes, fontsize=fontsize)
ax.text(0.45, 0.05, "Towards Inner Target", bbox=dict(color="white"), transform=ax.transAxes, fontsize=fontsize)
#ax.text(0.7, 0.3, r"$\Delta$" + "M = 0.3", color=colors[3], rotation=10, bbox=dict(color="white"), transform=ax.transAxes, fontsize=fontsize)

fig.tight_layout()
fig.show()
