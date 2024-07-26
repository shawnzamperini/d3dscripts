# Plot and compare the profiles from DIVIMP of a blobby scenario and a diffusive
# scenario.
import matplotlib.pyplot as plt
import oedge_plots


bpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-002a.nc"
dpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mrc-shifted-nodrift-2.nc"
bop = oedge_plots.OedgePlots(bpath)
dop = oedge_plots.OedgePlots(dpath)

ring = 92
s, nzb = bop.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
s, nzd = dop.along_ring(ring, "DDLIMS", charge="all", plot_it=False)

nzb = nzb / nzb.max()
nzd = nzd / nzd.max()

fig, ax = plt.subplots(figsize=(5,4))
ax.plot(s, nzb, lw=3, label="Blob-induced")
ax.plot(s, nzd, lw=3, label="Diffusive")
ax.legend()
ax.set_xlabel("Distance from inner target (m)", fontsize=14)
ax.set_ylabel("W Density (normalized)", fontsize=14)
fig.tight_layout()
fig.show()
