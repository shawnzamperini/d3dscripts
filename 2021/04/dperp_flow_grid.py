# Script to plot a contour plot of max density values for different combinations
# of dperp and flow values.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker


# Constants for the plotting.
ring = 17; ring_end = 52

# Values defining the grid locations.
flows = [0.00, 0.08, 0.17, 0.25, 0.33, 0.42, 0.50]
dperps = [0.10, 0.16, 0.27, 0.45, 0.74, 1.21, 2.00]

# First load in all the cases.
case_root = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/dperp_flow_grid/d3d-167277-inj-grid"
cases = {}
for flow_num in range(1, len(flows) + 1):
    for dperp_num in range(1, len(dperps) + 1):
        case_num = str(flow_num) + str(dperp_num)
        print("Loading: {}".format(case_num))
        case_path = case_root + case_num + ".nc"
        cases[case_num] = oedge_plots.OedgePlots(case_path)

# For each case we want the maximum W density value along the ring (excluding
# the value at the injection location).
nzs = []; nz_maxes = {}
fig, ax = plt.subplots()
for case_num, op in cases.items():
    print("Along ring: {}".format(case_num))
    s, nz = op.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    nzs.append(nz)
    keep = s < ring_end
    nz_maxes[case_num] = nz[keep].max() / nz[keep].mean()

    # Optional plot.
    ax.plot(s, nz, label=case_num)
ax.set_xlabel("S (m)")
ax.set_ylabel("Density")
fig.tight_layout()
fig.show()

# Create mesh grid of the X, Y and Z values for plots.
X, Y = np.meshgrid(dperps, flows)
Z = np.reshape(list(nz_maxes.values()), (len(flows), len(dperps)))

# x, y pairs for each point location.
x = np.tile(dperps, len(flows))
y = np.repeat(flows, len(dperps))

# For the xticks.
xticklabels = [0.1, 0.5, 1.0, 2.0]

# Plotting constants, parameters.
plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)
fontsize = 16

fig, ax1 = plt.subplots(figsize=(8, 8*(5/8)))

cont = ax1.contourf(X, Y, Z, cmap="magma", levels=np.arange(0, 11, 1), extend="max")
cbar = fig.colorbar(cont)
#ax1.scatter(x, y, color="k")
ax1.set_xlabel(r"$\mathdefault{D}_{\perp}\ \mathdefault{(m^2/s)}$", fontsize=fontsize)
ax1.set_ylabel("Imposed Mach Flow", fontsize=fontsize)
cbar.ax.set_ylabel(r"$\mathdefault{n_z^{max}/n_z^{average}}$", fontsize=fontsize)
cbar.ax.tick_params(axis="both", which="both", labelsize=14)
ax1.set_xscale("log")
ax1.set_xticks(xticklabels)
ax1.set_xticklabels(xticklabels)
ax1.tick_params(axis="both", which="both", labelsize=14)

fig.tight_layout()
fig.show()
