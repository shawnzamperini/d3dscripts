# This script will plot the density and the criteria for accumulation.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Input parameters.
charge      = 15
#poly_order  = 11

path1 = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
path2 = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
op1 = oedge_plots.OedgePlots(path1)
op2 = oedge_plots.OedgePlots(path2)

ops = [op1, op2]
ss = []; vzs_raw = []; nzs_raw = []; gzs_raw = []; vzs_fit = []; nzs_fit = []
gzs_fit = []; dnzs_fit = []; crit3s = []
for i in range(0, len(ops)):

    # Case specific parameters.
    if ops[i] == op1:
        #stag_region = [15, 50]; play = 5; ring = 17; root_num = 0; poly_order = 7
        play = 0.5; ring = 17; window = 15; power = 5; fit_region = [15, 50]; zero = 42.60
    elif ops[i] == op2:
        play = 1; ring = 17; window = 15; power = 5; fit_region = [15, 50]; zero = 25.98+0.5

    # Savgol filters to density and velocity.
    s, vz = ops[i].along_ring(ring, "VELAVG", charge=charge, plot_it=False)
    s, nz = ops[i].along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
    keep = np.logical_and(s >= fit_region[0], s <= fit_region[1])
    s  = s[keep]
    vz = vz[keep]
    nz = nz[keep]
    vz_fit = savgol_filter(vz, window, power)
    nz_fit = savgol_filter(nz, window, power)
    gz_fit = nz_fit * vz_fit

    # Derivatives.
    dvz = np.gradient(vz_fit, s)
    dnz = np.gradient(nz_fit, s)
    dgz = np.gradient(gz_fit, s)
    dvz2 = np.gradient(dvz, s)
    dnz2 = np.gradient(dnz, s)
    dgz2 = np.gradient(dgz, s)

    # Criteria 3.
    crit3 = dvz2 / vz - dgz2 / gz_fit
    crit3 = savgol_filter(crit3, window, power)

    # Replace data near the sigularity with nans.
    crit3[np.logical_and(s>=zero-play, s<=zero+play)] = np.nan

    # Print the average values right outside the singularity.
    crit3_left = np.where(np.isnan(crit3))[0][0] - 1
    crit3_right = np.where(np.isnan(crit3))[0][-1] + 1
    crit3_avg = np.mean((crit3[crit3_left], crit3[crit3_right]))
    print("Average Criteria 3 = {:.3f}".format(crit3_avg))

    crit4 = savgol_filter(dgz2 / dvz2, window, power)
    fig, ax = plt.subplots()
    ax.plot(s, nz_fit)
    ax.plot(s, crit4)
    fig.tight_layout()
    fig.show()

    # Append to lists, normalized to the data.
    ss.append(s)
    vzs_raw.append(vz/np.abs(vz).max())
    nzs_raw.append(nz/np.abs(nz).max())
    vzs_fit.append(vz_fit/np.abs(vz).max())
    nzs_fit.append(nz_fit/np.abs(nz).max())
    gzs_fit.append(gz_fit/np.abs(gz_fit).max())
    dnzs_fit.append(dnz/np.abs(dnz).max())
    crit3s.append(crit3/np.nanmax(np.abs(crit3)))


# Distance to separatrix at OMP.
mid_dist = ops[0].nc.variables["MIDIST"][1][ring]
mid_str = "R-" + r"$\mathdefault{R_{sep}}$" + " = {:.2f} cm".format(mid_dist*100)

# Plotting.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5
combine = False

# The plots.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 8))

ax1.axhline(0, color="k")
ax1.plot(ss[0], dnzs_fit[0], color=colors[2], label=r"$\mathdefault{dn_Z}$/ds", lw=lw)
ax1.plot(ss[0], vzs_raw[0], color=colors[3], alpha=0.4, lw=lw)
ax1.plot(ss[0], vzs_fit[0], color=colors[3], label=r"$\mathdefault{v_Z}$", lw=lw)
ax1.plot(ss[0], nzs_raw[0], color=colors[1], alpha=0.4, lw=lw)
ax1.plot(ss[0], nzs_fit[0], color=colors[1], label=r"$\mathdefault{n_Z}$", lw=lw)
ax1.plot(ss[0], crit3s[0],  color=colors[4], label="Criteria", lw=lw)
ax1.legend(fontsize=12)
ax1.grid(zorder=1)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_ylabel("Normalized values", fontsize=fontsize)
#ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.text(0.05, 0.1, mid_str, transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.text(0.05, 0.95, "a) W15+ no additional flow", transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.set_ylim([-1, 1.2])
ax1.set_xlim([38, 46])

ax2.axhline(0, color="k")
ax2.plot(ss[1], dnzs_fit[1], color=colors[2], label=r"$\mathdefault{dn_Z}$/ds", lw=lw)
ax2.plot(ss[1], vzs_raw[1], color=colors[3], alpha=0.4, lw=lw)
ax2.plot(ss[1], vzs_fit[1], color=colors[3], label=r"$\mathdefault{v_Z}$", lw=lw)
ax2.plot(ss[1], nzs_raw[1], color=colors[1], alpha=0.4, lw=lw)
ax2.plot(ss[1], nzs_fit[1], color=colors[1], label=r"$\mathdefault{n_Z}$", lw=lw)
ax2.plot(ss[1], crit3s[1],  color=colors[4], label="Criteria", lw=lw)
ax2.grid(zorder=1)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.set_ylabel("Normalized values", fontsize=fontsize)
ax2.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax2.text(0.05, 0.95, "b) W15+ M = 0.4 additional flow", transform=ax2.transAxes, fontsize=fontsize, bbox=dict(color="white"))
#ax2.set_ylim([-0.6, 1.2])
ax2.set_ylim([-1.0, 1.2])

fig.tight_layout()
fig.show()
