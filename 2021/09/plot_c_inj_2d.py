import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge/")
import oedge_plots


ncpath = "/Users/zamperini/Documents/d3d_work/184267/d3d-184267-inj-005.nc"
op = oedge_plots.OedgePlots(ncpath)

op.plot_contour_polygon("DDLIMS", cmap="magma", cbar_label="C13 Density (m-3)",
  charge="all", scaling=op.absfac, normtype="log", vmin=5e13, vmax=5e15)
