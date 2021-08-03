# Plots of impurity density along a near-SOL field line.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Select series with or without drifts.
shot = 167247
case = "031"  # No drifts
#case = "032"  # Unfavorable ExB drifts 60%
#case = "033"  # Favorable ExB drifts 60%
ring = 30; inj_start = 58; xy = (62, 0.3)
force_ylim = [-3e-16, 3e-16]

shot = 167277
#case = "002"
case = "006"
ring = 17; inj_start = 65; xy = (69, 0.3)
force_ylim = [-4.5e-16, 4.5e-16]

# Paths to each netcdf file.
# m00 = No additional flow.
# m01 = Additional M = -0.1 flow
# etc...
unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}a.nc".format(167247, 167247, "034")
m00_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}.nc".format(shot, shot, case)  # No letter.
m01_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}a.nc".format(shot, shot, case)
m02_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}b.nc".format(shot, shot, case)
m03_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}c.nc".format(shot, shot, case)
m04_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}d.nc".format(shot, shot, case)
m05_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}e.nc".format(shot, shot, case)
mc00 = oedge_plots.OedgePlots(m00_path)
mc01 = oedge_plots.OedgePlots(m01_path)
mc02 = oedge_plots.OedgePlots(m02_path)
mc03 = oedge_plots.OedgePlots(m03_path)
mc04 = oedge_plots.OedgePlots(m04_path)
mc05 = oedge_plots.OedgePlots(m05_path)
unf  = oedge_plots.OedgePlots(unf_path)

# Load along ring total impurity density data.
mc00_along = list(mc00.along_ring(ring, "DDLIMS", charge="all", plot_it=False))
mc01_along = list(mc01.along_ring(ring, "DDLIMS", charge="all", plot_it=False))
mc02_along = list(mc02.along_ring(ring, "DDLIMS", charge="all", plot_it=False))
mc03_along = list(mc03.along_ring(ring, "DDLIMS", charge="all", plot_it=False))
mc04_along = list(mc04.along_ring(ring, "DDLIMS", charge="all", plot_it=False))
mc05_along = list(mc05.along_ring(ring, "DDLIMS", charge="all", plot_it=False))
unf_along  = list(unf.along_ring(24,    "DDLIMS", charge="all", plot_it=False))

# Get maximum impurity density value.
max_w = 0
for mc in [mc00_along, mc01_along, mc02_along, mc03_along, mc04_along, mc05_along]:
    if mc[1].data.max() > max_w:
        max_w = mc[1].data.max()

# Normalize.
for mc in [mc00_along, mc01_along, mc02_along, mc03_along, mc04_along, mc05_along]:
    mc[1] = mc[1] / max_w

# Normalize the unf case too.
unf_along[1] = unf_along[1] / np.max(unf_along[1])

# Get net force data, dropping s = 0 points.
charge = 8
vz_mult = 0.0
force = "fnet"
mc00_force = list(mc00.along_ring(ring, force, charge=charge, vz_mult=vz_mult, plot_it=False))
mc01_force = list(mc01.along_ring(ring, force, charge=charge, vz_mult=vz_mult, plot_it=False))
mc02_force = list(mc02.along_ring(ring, force, charge=charge, vz_mult=vz_mult, plot_it=False))
mc03_force = list(mc03.along_ring(ring, force, charge=charge, vz_mult=vz_mult, plot_it=False))
mc04_force = list(mc04.along_ring(ring, force, charge=charge, vz_mult=vz_mult, plot_it=False))
mc05_force = list(mc05.along_ring(ring, force, charge=charge, vz_mult=vz_mult, plot_it=False))
unf_force = list(unf.along_ring(24, force, charge=charge, vz_mult=vz_mult, plot_it=False))
mask = mc00_force[1] != 0.0
unf_mask = unf_force[1] != 0.0
mc00_force[0] = mc00_force[0][mask]; mc00_force[1] = mc00_force[1][mask]
mc01_force[0] = mc01_force[0][mask]; mc01_force[1] = mc01_force[1][mask]
mc02_force[0] = mc02_force[0][mask]; mc02_force[1] = mc02_force[1][mask]
mc03_force[0] = mc03_force[0][mask]; mc03_force[1] = mc03_force[1][mask]
mc04_force[0] = mc04_force[0][mask]; mc04_force[1] = mc04_force[1][mask]
mc05_force[0] = mc05_force[0][mask]; mc05_force[1] = mc05_force[1][mask]
unf_force[0] = unf_force[0][unf_mask]; unf_force[1] = unf_force[1][unf_mask]

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

stag_s = get_stag_s(unf, 24)
unf_x, unf_y = get_closest_imp(unf_along[0], unf_along[1], stag_s)
stag_s = get_stag_s(mc00, ring)
mc00_x, mc00_y = get_closest_imp(mc00_along[0], mc00_along[1], stag_s)
stag_s = get_stag_s(mc01, ring)
mc01_x, mc01_y = get_closest_imp(mc01_along[0], mc01_along[1], stag_s)
stag_s = get_stag_s(mc02, ring)
mc02_x, mc02_y = get_closest_imp(mc02_along[0], mc02_along[1], stag_s)
stag_s = get_stag_s(mc03, ring)
mc03_x, mc03_y = get_closest_imp(mc03_along[0], mc03_along[1], stag_s)
stag_s = get_stag_s(mc04, ring)
mc04_x, mc04_y = get_closest_imp(mc04_along[0], mc04_along[1], stag_s)
stag_s = get_stag_s(mc05, ring)
mc05_x, mc05_y = get_closest_imp(mc05_along[0], mc05_along[1], stag_s)

# Distance to separatrix at OMP.
mid_dist = mc00.nc.variables["MIDIST"][1][ring]
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(mid_dist*100)

# Plotting.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 6))
fontsize = 14
lw = 5
combine = False
stags = True
s = 100
zorder = 1
linewidths = 1.5

if combine:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 8))
else:
    fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.annotate("Unfavorable-"+r"$\mathdefault{B_T}$", (43.5, 1.0), xytext=(25, 1.15), color="blueviolet", bbox=dict(color="white"), arrowprops=dict(color="blueviolet", facecolor="blueviolet", arrowstyle="-"), fontsize=fontsize)
ax1.plot(unf_along[0],  unf_along[1], color="blueviolet", lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(unf_x, unf_y, s=s, edgecolor="k", color="blueviolet", zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.plot(mc00_along[0], mc00_along[1], label=r"$\Delta$" + "M = 0.0", color=colors[0], lw=lw, zorder=zorder)
#ax1.plot(mc00_along[0], mc00_along[1], label=r"$\Delta$" + "M = 0.0", color="lightcoral", lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(mc00_x, mc00_y, s=s, edgecolor="k", color=colors[0], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.plot(mc01_along[0], mc01_along[1], label=r"$\Delta$" + "M = 0.1", color=colors[1], lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(mc01_x, mc01_y, s=s, edgecolor="k", color=colors[1], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.plot(mc02_along[0], mc02_along[1], label=r"$\Delta$" + "M = 0.2", color=colors[2], lw=lw, zorder=zorder)
#ax1.plot(mc02_along[0], mc02_along[1], label=r"$\Delta$" + "M = 0.2", color="indianred", lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(mc02_x, mc02_y, s=s, edgecolor="k", color=colors[2], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.plot(mc03_along[0], mc03_along[1], label=r"$\Delta$" + "M = 0.3", color=colors[3], lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(mc03_x, mc03_y, s=s, edgecolor="k", color=colors[3], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.plot(mc04_along[0], mc04_along[1], label=r"$\Delta$" + "M = 0.4", color=colors[4], lw=lw, zorder=zorder)
#ax1.plot(mc04_along[0], mc04_along[1], label=r"$\Delta$" + "M = 0.4", color="maroon", lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(mc04_x, mc04_y, s=s, edgecolor="k", color=colors[4], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.plot(mc05_along[0], mc05_along[1], label=r"$\Delta$" + "M = 0.5", color=colors[5], lw=lw, zorder=zorder)
zorder += 1
if stags:
    ax1.scatter(mc05_x, mc05_y, s=s, edgecolor="k", color=colors[5], zorder=zorder, linewidths=linewidths)
    zorder += 1
ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.set_ylabel("Tungsten density (normalized)", fontsize=fontsize)
ax1.legend(fontsize=12, framealpha=1.0, loc=(0.025, 0.31))
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.grid()
ax1.set_ylim([0.0, 1.2])
ax1.axvline(inj_start, color="k", linestyle="--", lw=lw)
ax1.text(xy[0], xy[1], "W Injection Location", rotation=90, fontsize=fontsize, bbox=dict(color="white"))
#ax1.text(0.05, 0.95, mid_str, transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.text(0.05, 0.87, mid_str, transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
#ax1.text(0.025, 0.95, "a)", transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))

# Show the current figure and start a new one.
if not combine:
    fig.tight_layout()
    fig.show()
    fig, ax2 = plt.subplots(figsize=(5, 4))

ax2.axhline(0.0, color="k", lw=3, linestyle="--")
ax2.plot(unf_force[0], unf_force[1], label="Unfavorable", color="m", lw=lw)
ax2.plot(mc00_force[0], mc00_force[1], label="M = 0.0", color=colors[0], lw=lw)
ax2.plot(mc01_force[0], mc01_force[1], label="M = 0.1", color=colors[1], lw=lw)
ax2.plot(mc02_force[0], mc02_force[1], label="M = 0.2", color=colors[2], lw=lw)
ax2.plot(mc03_force[0], mc03_force[1], label="M = 0.3", color=colors[3], lw=lw)
ax2.plot(mc04_force[0], mc04_force[1], label="M = 0.4", color=colors[4], lw=lw)
ax2.plot(mc05_force[0], mc05_force[1], label="M = 0.5", color=colors[5], lw=lw)
ax2.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax2.set_ylabel("Net Force W{}+ (N)".format(charge), fontsize=fontsize)
ax2.legend(fontsize=12, framealpha=1.0, loc="upper right", ncol=2)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.grid()
ax2.set_ylim(force_ylim)
ax2.text(0.33, 0.67, r"$\mathdefault{F_{net}}$", transform=ax2.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax2.arrow(0.26, 0.62, 0.2, 0.0, transform=ax2.transAxes, width=0.01, color="k", zorder=50)
ax2.text(0.65, 0.17, r"$\mathdefault{F_{net}}$", transform=ax2.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax2.arrow(0.8, 0.25, -0.2, 0.0, transform=ax2.transAxes, width=0.01, color="k", zorder=50)
ax2.text(0.025, 0.95, "b)", transform=ax2.transAxes, fontsize=fontsize, bbox=dict(color="white"))

fig.tight_layout()
fig.show()

# Find where the density peaks. Assumes ring 17 on 167277.
s_peaks = []
for mc in [mc00_along, mc01_along, mc02_along, mc03_along, mc04_along, mc05_along]:
    keep = np.logical_and(mc[0] > 10, mc[0] < 50)
    max_idx = mc[1][keep].argmax()
    s_peaks.append(mc[0][keep][max_idx])
stag_points = np.array([42.6, 42.0, 38.7, 29.9, 26.1, 22.6])
zero_force = np.array([42.5, 42.2, 42.07, 38.23, 33.44, 28.62])
machs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

fig, ax3 = plt.subplots(figsize=(5, 4))

ax3.scatter(machs, s_peaks-stag_points)
#ax3.scatter(machs, s_peaks-zero_force)
ax3.set_xlabel("Imposed Mach flow")
ax3.set_ylabel("Distance from peak density")

fig.tight_layout()
fig.show()
