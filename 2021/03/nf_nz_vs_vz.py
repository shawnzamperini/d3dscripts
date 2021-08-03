# The script is a plot of the location of the stagnation point vs the location
# of the peak impurity density.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


charge = 9
ring   = 17

ops = []
for case in ["", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]:
    path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006{}.nc".format(case)
    ops.append(oedge_plots.OedgePlots(path))

machs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.05, 0.15, 0.25, 0.35, 0.45]

# Values for W9+ on ring 17 for each case in machs.
nz_peaks = [49.4, 49.2, 49.03, 49.06]
vz_stags = [41.1, 39.0, 31.22, 26.00]

nz_peaks = np.array([]); vz_stags = np.array([])
for op in ops:
    s, nz = op.along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
    s, vz = op.along_ring(ring, "VELAVG", charge=charge, plot_it=False)
    keep = np.logical_and(s > 10, s < 50)
    max_idx = nz[keep].argmax()
    stag_idx = np.abs(vz[keep]).argmin()
    nz_peaks = np.append(nz_peaks, s[keep][max_idx])
    vz_stags = np.append(vz_stags, s[keep][stag_idx])

fig, ax1 = plt.subplots()

ax1.axhline(0, color="k", linestyle="--")
ax1.scatter(machs, vz_stags-nz_peaks)

fig.tight_layout()
fig.show()


# The following snippet plots the impurity velocity and density so that we can
# get values for nz_peaks and vz_stags.
i = 3
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
s, nz = ops[i].along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
s, vz = ops[i].along_ring(ring, "VELAVG", charge=charge, plot_it=False)
ax1.axhline(0, color="k")
ax2.axhline(0, color="k")
ax1.plot(s, nz)
ax2.plot(s, vz)
ax1.set_ylabel("DDLIMS")
ax2.set_ylabel("VELAVG")
fig.tight_layout()
fig.show()
