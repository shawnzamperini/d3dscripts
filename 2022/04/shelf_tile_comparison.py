# Script to show some metrics comparing graphite to SiC shelf tiles.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter


# Load runs.
print("Loading runs...")
grap_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-grap-shelf.nc"
sic_c_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-sic-c-shelf.nc"
sic_si_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-sic-si-shelf.nc"
grap = oedge_plots.OedgePlots(grap_path)
sic_c = oedge_plots.OedgePlots(sic_c_path)
sic_si = oedge_plots.OedgePlots(sic_si_path)

# Load along ring profiles of the respective impurity density.
ring = 20
sgr, nzgr = grap.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
ssc, nzsc = sic_c.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
sss, nzss = sic_si.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
nzgrs = savgol_filter(nzgr, 25, 2)
nzscs = savgol_filter(nzsc, 25, 2)
nzsss = savgol_filter(nzss, 25, 2)

# See what these tiles' contribution to Zeff is.
sgr, ne = grap.along_ring(ring, "KNBS", plot_it=False)
grap_zeff = ne.copy()  # Approximating here that nD ~ ne
sic_zeff = ne.copy()
#grap_tot_dens = ne
#sic_tot_dens = ne

# The carbon ions first for both cases.
for i in range(0, 6):
    sgr, nzgr_tmp = grap.along_ring(ring, "DDLIMS", charge=i+1, plot_it=False)
    ssc, nzsc_tmp = sic_c.along_ring(ring, "DDLIMS", charge=i+1, plot_it=False)
    #grap_tot_dens += nzgr_tmp
    #sic_tot_dens += nzsc_tmp
    grap_zeff += nzgr_tmp * (i + 1)**2
    sic_zeff += nzsc_tmp * (i + 1)**2
grap_zeff = grap_zeff / ne

# Do aagin for the Si ions in the SiC case.
for i in range(0, 14):
    sss, nzss_tmp = sic_si.along_ring(ring, "DDLIMS", charge=i+1, plot_it=False)
    #sic_tot_dens += nzss_tmp
    sic_zeff += nzss_tmp * (i + 1)**2
sic_zeff = sic_zeff / ne

grap_zeffs = savgol_filter(grap_zeff, 25, 2)
sic_zeffs = savgol_filter(sic_zeff, 25, 2)

rmrsomp = grap.nc.variables["MIDIST"][1][ring]*1000

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

alpha = 1.0
ax1.plot(sgr, nzgr, color="k", label="Graphite", alpha=alpha)
ax1.plot(ssc, nzsc, color="k", label="SiC (C)", linestyle="--", alpha=alpha)
ax1.plot(sss, nzss, color="r", label="SiC (Si)", linestyle="--", alpha=alpha)
#ax1.plot(sgr, nzgrs, label="Graphite", color="k")
#ax1.plot(ssc, nzscs, label="SiC (C)", color="r")
#ax1.plot(sss, nzsss, label="SiC (Si)", color="r", linestyle="--")
ax1.legend()
ax1.set_ylim([1e14, 1.2e17])
ax1.set_yscale("log")
ax1.set_xlabel("Distance from inner target (m)")
ax1.set_ylabel("Impurity density (m-3)")

ax2.axhline(1.0, color="k", linestyle="--")
ax2.plot(sgr, grap_zeff, label="Graphite", color="k",alpha=alpha)
ax2.plot(ssc, sic_zeff, label="SiC", color="r", alpha=alpha)
#ax2.plot(sgr, grap_zeffs, label="Graphite", color="k")
#ax2.plot(ssc, sic_zeffs, label="SiC", color="r")
ax2.set_xlabel("Distance from inner target (m)")
ax2.set_ylabel("Zeff")
ax2.legend()
ax2.set_ylim([0.95, 1.30])

fig.suptitle("#167196, R-Rsep OMP = {:.2f} mm".format(rmrsomp))
fig.tight_layout()
fig.show()
