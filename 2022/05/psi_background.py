# Script to show a plot fo the background plasma.
import oedge_plots
import matplotlib.pyplot as plt



ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/184267/d3d-184267-inj-015.nc"
plt.rcParams["font.family"] = "Century Gothic"

op = oedge_plots.OedgePlots(ncpath)

# Te plot.
fig = op.plot_contour_polygon("KTEBS", cmap="inferno",
    cbar_label=r"$\mathdefault{T_e\ (eV)}$", normtype="log", figsize=(5, 6),
    vmin=9, vmax=100, show_grid=False)
fig.axes[0].set_axis_off()
fig.show()

# 13C plot.
fig = op.plot_contour_polygon("DDLIMS", cmap="nipy_spectral", charge="all",
    cbar_label=r"$\mathdefault{^{13}C\ Density\ (m^{-3})}$", normtype="log",
    figsize=(5, 6), show_grid=False, scaling=op.absfac, vmin=1e14, vmax=1e16)
fig.axes[0].set_axis_off()
fig.show()
