# Return the surfaced averaged values of tungsten concentration from our
# blob DIVIMP sim.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt


ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div6.nc"
op = oedge_plots.OedgePlots(ncpath)
#aminor = 0.588
aminor = (2.260 - 1.723)
r0 = 1.723

core_rings = np.arange(5, 19)
cws = []
rhos = []
for ring in core_rings:
    s, nz = op.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    s, ne = op.along_ring(ring, "KNBS", plot_it=False)
    cws.append((nz/ne).mean())

    # Calculate rho.
    s, r = op.along_ring(ring, "RS", plot_it=False)
    romp = r.max()
    rho = (romp - r0) / aminor
    rhos.append(rho)

# Load the SXR data.
path = "/Users/zamperini/Documents/d3d_work/files/imp_analysis_167196.npz"
imps = np.load(path)
rho = imps["rho"]
cw = imps["cw"]
cw_avg = np.nanmean(cw, axis=0)
cw_std = np.nanstd(cw, axis=0)

fig, ax1 = plt.subplots()

ax1.plot(rho, cw_avg, color="r")
ax1.fill_between(rho, cw_avg-cw_std, cw_avg+cw_std, color="r", alpha=0.3)
ax1.scatter(rhos, cws, s=25, marker="*", color="r")

ax1.set_xlim([0.7, 1])
ax1.set_ylim([0, 0.005])
ax1.set_xlabel("Rho")
ax1.set_ylabel("W Concentration")
fig.tight_layout()
fig.show()
