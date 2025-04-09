import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


# Load files
root = "/fusion/projects/codes/oedge/zamperinis/results"
ncpath = "{}/d3d-w-wall-param-scan-v6-w-all-1.nc".format(root)
op_pr = oedge_plots.OedgePlots(ncpath)
ncpath = "{}/d3d-w-wall-param-scan-v6-w-all-1-noprompt.nc".format(root)
op_np = oedge_plots.OedgePlots(ncpath)

# Concentrations
ring = 15
s, nz_pr = op_pr.along_ring(ring, "DDLIMS", charge="all", scaling=op_pr.absfac, plot_it=False)
s, nz_np = op_np.along_ring(ring, "DDLIMS", charge="all", scaling=op_np.absfac, plot_it=False)
s, ne_s = op_pr.along_ring(ring, "KNBS", plot_it=False)
cz_pr = nz_pr / (nz_pr + ne_s)
cz_np = nz_np / (nz_np + ne_s)

# Plot
fig, ax = plt.subplots(figsize=(5, 4))
#ax.plot(s, cz_np, lw=3, color="k")
#ax.plot(s, cz_np, label="OFF", lw=2, color="tab:red")
ax.plot(s, cz_pr, lw=3, color="k")
ax.plot(s, cz_pr, label="ON", lw=2, color="tab:purple")
ax.set_yscale("log")
ax.set_ylim([1e-3, 1e-1])
ax.set_xlabel("Distance from inner target (m)", fontsize=16)
ax.set_ylabel("W Concentration", fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=14)
#ax.legend()
ax.grid(alpha=0.5, which="both")
fig.tight_layout()
fig.show()
