#
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np

unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/d3d-{}-inj-{}a.nc".format(167247, 167247, "034")
op = oedge_plots.OedgePlots(unf_path)
op.plot_contour_polygon("DDLIMS", charge="all", cmap="magma", vmin=0.001, vmax=0.003)
