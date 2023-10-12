# Quick script to plot the 2D W density from a blobby transport model.
import oedge_plots

# 020b - parallel transport OFF
# 021  - parallel transport ON
#blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-020b-nocore.nc"
#blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-021.nc"
blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-fluc-002.nc"
blob = oedge_plots.OedgePlots(blob_nc_path)

fig = blob.plot_contour_polygon("DDLIMS", scaling=blob.absfac, cmap="nipy_spectral", vmin=1e11, vmax=1e15,
                                cbar_label=r"W Density $\mathdefault{(m^{-3})}$", show_mr=True, charge="all",
                                normtype="log", show_cp=True, ptip=2.282)

ax = fig.axes[0]
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)