import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


# Load run into OedgePlots object.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/oedge_tutorial/d3d-167196-divimp-csput-v1.nc"
op = oedge_plots.OedgePlots(ncpath)

# 2D plot of the W density (all charge states).  DDLIMS is the name of impurity density array in netCDF file. It
# must be scaled by op.absfac to go from normalized to physical units.
op.plot_contour_polygon("DDLIMS", charge="all", normtype="log", cmap="inferno", lut=7, vmin=1e9, vmax=5e15,
                        scaling=op.absfac, cbar_label="W Density (m-3)")

# 2D plot of just the W15+ density.
op.plot_contour_polygon("DDLIMS", charge=5, normtype="log", cmap="inferno", lut=7, vmin=1e9, vmax=5e15,
                        scaling=op.absfac, cbar_label="W5+ Density (m-3)")

# 2D plot of the W concentration. We make use of the own_data parameter to pass in manipulated data. We can mask the
# data so it doesn't plot anything where it equals zero.
wdens = op.read_data_2d("DDLIMS", charge="all", scaling=op.absfac)
ne = op.read_data_2d("KNBS")
wconc = np.ma.masked_where(wdens==0, wdens / ne)
op.plot_contour_polygon(own_data=wconc, normtype="log", lut=5, cbar_label="W Concentration", cmap="cool",
                        vmin=1e-8, vmax=1e-3)

# Plot of the W density along the separatrix ring.
s, nw = op.along_ring(op.irsep, "DDLIMS", charge="all", plot_it=False, scaling=op.absfac)
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(s, nw, color="k", lw=3)
ax.plot(s, nw, color="tab:red", lw=2)
ax.set_yscale("log")
ax.set_xlabel("Distance from inner target (m)", fontsize=12)
ax.set_ylabel("W Density (m-3)", fontsize=12)
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.show()

# Plot of the W density at the outboard midplane. Zeros are masked.
w_dict = op.along_line(2.18, 2.30, op.z0, op.z0, "DDLIMS", charge="all", scaling=op.absfac)
mask = np.array(w_dict["DDLIMS"]) > 0
x = np.array(w_dict["psin"])[mask]
y = np.array(w_dict["DDLIMS"])[mask]
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(x, y, color="k", lw=3)
ax.plot(x, y, color="tab:red", lw=2)
ax.set_yscale("log")
ax.set_xlabel("Psin", fontsize=12)
ax.set_ylabel("W Density (m-3)", fontsize=12)
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.show()