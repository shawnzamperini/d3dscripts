# Script to plot near-SOL C13 profiles from DIVIMP.
import oedge_plots
import matplotlib.pyplot as plt


include_lsn = True
lsn_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d-167196-modE-cd4puff-lowerosp-d03.nc"
unf_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184527/d3d-184527-inj-015.nc"
fav_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-inj-015.nc"
unf = oedge_plots.OedgePlots(unf_path)
fav = oedge_plots.OedgePlots(fav_path)
lsn = oedge_plots.OedgePlots(lsn_path)

ring = 15
lsn_ring = 24
sunf, nzunf = unf.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
sfav, nzfav = fav.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
slsn, nzlsn = lsn.along_ring(lsn_ring, "DDLIMS", charge="all", plot_it=False)
nzlsn = nzlsn * unf.absfac
slsn = slsn[::-1]
nzlsn = nzlsn[::-1]
print("MET: R-Rsep OMP = {:.2f}".format(unf.nc.variables["MIDIST"][1][ring]*1000))
print("LSN: R-Rsep OMP = {:.2f}".format(lsn.nc.variables["MIDIST"][1][lsn_ring]*1000))
# If we want to include the LSN results then we gotta normalize the s coordinate
# first.
if include_lsn:
    sfav = sfav / sfav.max()
    sunf = sunf / sunf.max()
    slsn = slsn / slsn.max()

fontsize = 14
fig, ax = plt.subplots(figsize=(5, 4))
if include_lsn:
    #ax.plot(slsn, nzlsn, lw=3, color="k", linestyle="--", label="LSN")
    ax.plot(slsn, nzlsn, lw=3, color="k", linestyle="--", label="LSN (Open)")
    ax.set_xlabel("Parallel distance from\nouter target (normalized)", fontsize=fontsize)
    ax.set_yscale("log")
    ax.grid()
    ax.set_ylim([1e14, 4e17])
else:
    ax.set_xlabel("Distance from outer target (m)", fontsize=fontsize)
#ax.plot(sunf, nzunf, color="tab:purple", lw=3, label="Unfavorable")
#ax.plot(sfav, nzfav, color="tab:red", lw=3, label="Favorable")
ax.plot(sfav, nzfav, color="tab:red", lw=3, label="USN (Closed)")
ax.legend(fontsize=fontsize-2)
ax.tick_params(which="both", labelsize=14)

ax.set_ylabel(r"$\mathdefault{^{13}C\ Density\ (m^{-3})}$", fontsize=fontsize)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()

ring = 36
sunf, nzunf = unf.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
sfav, nzfav = fav.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
print("R-Rsep OMP = {:.2f}".format(unf.nc.variables["MIDIST"][1][ring]*1000))

fontsize = 14
fig, ax = plt.subplots()
ax.plot(sunf, nzunf, color="tab:purple", lw=3, label="Unfavorable")
ax.plot(sfav, nzfav, color="tab:red", lw=3, label="Favorable")
ax.legend(fontsize=fontsize)
ax.set_xlabel("Distance from outer target (m)", fontsize=fontsize)
ax.set_ylabel(r"$\mathdefault{^{13}C\ Density\ (m^{-3})}$", fontsize=fontsize)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
