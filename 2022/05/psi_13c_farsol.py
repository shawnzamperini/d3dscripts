# Script to plot the far-SOL distribution of 13C.
import oedge_plots
import matplotlib.pyplot as plt


fav_path = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-inj-015.nc"
fav = oedge_plots.OedgePlots(fav_path)


ring = 36
sfav, nzfav = fav.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
print("R-Rsep OMP = {:.2f}".format(fav.nc.variables["MIDIST"][1][ring]*1000))

fontsize = 14
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(sfav, nzfav, color="tab:red", lw=3, label="Favorable")
#ax.legend(fontsize=fontsize)
ax.set_xlabel("Distance from outer target (m)", fontsize=fontsize)
ax.set_ylabel(r"$\mathdefault{^{13}C\ Density\ (m^{-3})}$", fontsize=fontsize)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
