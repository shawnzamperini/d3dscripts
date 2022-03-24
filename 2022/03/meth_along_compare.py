# Meh, simple script to just look at the along ring porfiles of the C13 density.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt


# Case names.
nc267_1 = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-inj-012.nc"
nc267_2 = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-inj-013.nc"
nc527_1 = "/Users/zamperini/Documents/d3d_work/divimp_files/184527/d3d-184527-inj-012.nc"
ring = 17
labels = ["F: No Flow", "F: M = 0.3", "U: No Flow"]
ncpaths = [nc267_1, nc267_2, nc527_1]

ops = []; ss = []; nzs = []
for i in range(len(ncpaths)):
    op = oedge_plots.OedgePlots(ncpaths[i])
    s, nz = op.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    ops.append(op)
    ss.append(s)
    nzs.append(nz)

fig, ax = plt.subplots()
for i in range(0, len(ops)):
    ax.plot(ss[i], nzs[i], label=labels[i])
ax.set_xlabel("Distance from outer target (m)")
ax.set_ylabel("C13 Density (m-3)")
ax.legend()
fig.tight_layout()
fig.show()
