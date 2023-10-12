import oedge_plots
import matplotlib.pyplot as plt


opd = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-018d.nc")
opg = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-018g.nc")
op_good = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010f.nc")

s, nwd = opd.along_ring(91, "DDLIMS", charge="all", plot_it=False)
s, nwg = opg.along_ring(91, "DDLIMS", charge="all", plot_it=False)
s, nw_good = op_good.along_ring(91, "DDLIMS", charge="all", plot_it=False)

fig, ax = plt.subplots()
ax.plot(s, nwd/nwd.max(), label="018d")
ax.plot(s, nwg/nwg.max(), label="018g")
ax.plot(s, nw_good/nw_good.max(), label="010f")
ax.legend()
ax.set_xlabel("s (m)")
ax.set_ylabel("nW norm")
fig.tight_layout()
fig.show()
