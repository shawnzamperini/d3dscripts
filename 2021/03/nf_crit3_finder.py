# This script will help find each runs criteria 3 value.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Input parameters.
charge      = 15

path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
op = oedge_plots.OedgePlots(path)

ops = [op]
ss = []; vzs_raw = []; nzs_raw = []; gzs_raw = []; vzs_fit = []; nzs_fit = []
gzs_fit = []; dnzs_fit = []; crit3s = []
for i in range(0, len(ops)):

    # Case specific parameters.
    if path == "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc":
        play = 0.0; ring = 17; window = 11; power = 2; fit_region = [35, 50]; zero = 42.60  # tuned
    elif path == "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006a.nc":
        play = 0.0; ring = 17; window = 11; power = 2; fit_region = [30, 46]; zero = 41.20
    elif path == "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006b.nc":
        play = 0.0; ring = 17; window = 21; power = 2; fit_region = [25, 50]; zero = 42.60
    elif path == "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006c.nc":
        play = 0.0; ring = 17; window = 31; power = 2; fit_region = [20, 50]; zero = 42.60
    elif path == "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc":
        #play = 1; ring = 17; window = 21; power = 3; fit_region = [20, 35]; zero = 25.98+0.5  # tuned
        play = 0; ring = 17; window = 5; power = 3; fit_region = [20, 35]; zero = 25.98+0.5  # tuned
    elif path == "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006e.nc":
        play = 0.5; ring = 17; window = 11; power = 2; fit_region = [35, 50]; zero = 42.60

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
    #crit3 = savgol_filter(crit3, window, power)

    # Replace data near the sigularity with nans.
    if play != 0.0:
        crit3[np.logical_and(s>=zero-play, s<=zero+play)] = np.nan

    # Print the average values right outside the singularity.
    if play != 0.0:
        crit3_left = np.where(np.isnan(crit3))[0][0] - 1
        crit3_right = np.where(np.isnan(crit3))[0][-1] + 1
        crit3_avg = np.mean((crit3[crit3_left], crit3[crit3_right]))
        print("Average Criteria 3 = {:.3f}  ({:.3f}  {:.3f})".format(crit3_avg, crit3[crit3_left], crit3[crit3_right]))
    else:

        # As defined crit 3 has to happen when dn/ds = 0.
        crit3_val = crit3[np.argmax(nz_fit)]
        print("Criteria 3 = {:.3f}".format(crit3_val))

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
fig, ax1 = plt.subplots(figsize=(5, 5))

ax1.axhline(0, color="k")
ax1.plot(ss[0], dnzs_fit[0], color=colors[2], label=r"$\mathdefault{dn_Z}$/ds", lw=lw)
ax1.plot(ss[0], vzs_raw[0], color=colors[3], alpha=0.4, lw=lw)
ax1.plot(ss[0], vzs_fit[0], color=colors[3], label=r"$\mathdefault{v_Z}$", lw=lw)
ax1.plot(ss[0], nzs_raw[0], color=colors[1], alpha=0.4, lw=lw)
ax1.plot(ss[0], nzs_fit[0], color=colors[1], label=r"$\mathdefault{n_Z}$", lw=lw)
ax1.plot(ss[0], crit3s[0],  color=colors[4], label="Criteria", lw=lw)
#ax1.plot(ss[0], gzs_fit[0],  color=colors[4], label="gamma", lw=lw)
ax1.legend(fontsize=12)
ax1.grid(zorder=1)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_ylabel("Normalized values", fontsize=fontsize)
#ax1.set_xlabel("Distance from inner target (m)", fontsize=fontsize)
ax1.text(0.05, 0.1, mid_str, transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.text(0.05, 0.95, "a) W15+ no additional flow", transform=ax1.transAxes, fontsize=fontsize, bbox=dict(color="white"))
ax1.set_ylim([-1, 1.2])
#ax1.set_xlim([38, 46])

fig.tight_layout()
fig.show()
