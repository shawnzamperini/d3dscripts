# This script will compare different values of Dperp on impurity accumulation.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
import matplotlib as mpl

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

#ring = 20; inj_start = 58; xy = (62, 0.3); xlim = [None, None]
ring = 30; inj_start = 58; xy = (62, 0.3); xlim = [None, None]
#ring = 18; inj_start = 82.5; xy = (86.5, 0.3); xlim = [20, 90]
#ring = 22; inj_start = 68.5; xy = (72.5, 0.3); xlim = [10, 80]
#ring = 30; inj_start = 68.5; xy = (72.5, 0.3); xlim = [10, 80]
#charge = 15
charge = "all"

dperps = [0.3, 0.4, 0.6, 1.0, 5.0]
dp03_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a2.nc"
dp06_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034g2.nc"
dp1_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034b2.nc"
dp5_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034c2.nc"
dp10_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034d2.nc"

dp03 = oedge_plots.OedgePlots(dp03_path)
dp06 = oedge_plots.OedgePlots(dp06_path)
dp1 = oedge_plots.OedgePlots(dp1_path)
dp5 = oedge_plots.OedgePlots(dp5_path)
dp10 = oedge_plots.OedgePlots(dp10_path)

def get_stag_s(op, ring):
    charges = np.arange(0, 31)
    num_knots = int(op.nc.variables["MAXNKS"][:])
    vzs = np.zeros((len(charges), num_knots))
    nzs = np.zeros((len(charges), num_knots))
    for j in range(0, len(charges)):
        charge = charges[j]
        s, vz = op.along_ring(ring, "VELAVG", charge=charge, plot_it=False, remove_zeros=False)
        s, nz = op.along_ring(ring, "DDLIMS", charge=charge, plot_it=False, remove_zeros=False)
        vzs[j] = vz
        nzs[j] = nz

    vzs_weighted = np.zeros(num_knots)
    for j in range(0, num_knots):
        if nzs[:, j].sum() != 0.0:
            vzs_weighted[j] = np.average(vzs[:, j], weights=nzs[:, j])

    keep = np.logical_and(vzs_weighted != 0.0, np.logical_and(s>=10, s<=55))
    vzs_weighted = vzs_weighted[keep]
    s = s[keep]
    stag_idx = np.abs(vzs_weighted).argmin()
    return s[stag_idx]

def get_closest_imp(x, y, s):
    idx = np.abs(x - s).argmin()
    return (x[idx], y[idx])

max_w = 0.0; all_imps = []
for op in [dp03, dp06, dp1, dp5, dp10]:
    s, imps = list(op.along_ring(ring, "DDLIMS", charge=charge, plot_it=False))
    if imps.data.max() > max_w:
        max_w = imps.data.max()
    all_imps.append(imps.data)

# Normalize.
for i in range(0, len(all_imps)):
    #all_imps[i] = all_imps[i] / max_w
    all_imps[i] = all_imps[i] / all_imps[i].max()
    #pass

ops = [dp03, dp06, dp1, dp5, dp10]
xs = []; ys = []
for i in range(0, len(ops)):
    stag_s = get_stag_s(ops[i], ring)
    x, y = get_closest_imp(s, all_imps[i], stag_s)
    xs.append(x)
    ys.append(y)

# Distance to separatrix at OMP.
mid_dist = dp03.nc.variables["MIDIST"][1][ring]
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(mid_dist*100)

# Find radial gradient of gamma at specified knot (should align with
# accumulation location).
knot = 87  # Can fine tune if needed.
charge = 15

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5
size = 100
zorder = 1
linewidths = 1.5

fig, ax1 = plt.subplots(figsize=(5, 4))

for i in range(0, len(all_imps)):
    label = r"$\mathdefault{D}_{\perp}$" + "= {:.2f}".format(dperps[i])
    ax1.plot(s, all_imps[i], label=label, color=colors[i], lw=lw, zorder=zorder)
    zorder += 1
    ax1.scatter(xs[i], ys[i], s=size, edgecolor="k", color=colors[i], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.set_ylabel("Tungsten density (normalized)", fontsize=fontsize)
ax1.grid(zorder=0)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.legend(fontsize=12)
ax1.axvline(inj_start, color="k", linestyle="--", lw=lw)
ax1.text(xy[0], xy[1], "W Injection Location", rotation=90, fontsize=fontsize, bbox=dict(color="white"))
ax1.set_ylim([0.0, 1.2])
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_xlim(xlim)
ax1.text(0.4, 0.92, mid_str, fontsize=fontsize, bbox=dict(color="white"), transform=ax1.transAxes)

fig.tight_layout()
fig.show()
