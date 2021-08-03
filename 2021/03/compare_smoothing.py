# Some plots to compare levels of smoothing in OSM.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
import matplotlib as mpl

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

op0_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-bkg-005_smooth0.nc"
op1_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-bkg-005_smooth1.nc"
op2_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-bkg-005_smooth2.nc"
op0 = oedge_plots.OedgePlots(op0_path)
op1 = oedge_plots.OedgePlots(op1_path)
op2 = oedge_plots.OedgePlots(op2_path)


# Compare Ti for a near-SOL ring.
fontsize = 14
s, ti0 = op0.along_ring(40, "KTIBS", plot_it=False)
s, ti1 = op1.along_ring(40, "KTIBS", plot_it=False)
s, ti2 = op2.along_ring(40, "KTIBS", plot_it=False)
fig, ax1 = plt.subplots()
ax1.plot(s, ti0, label="278=0")
ax1.plot(s, ti1, label="278=1")
ax1.plot(s, ti2, label="278=2")
ax1.legend(fontsize=12)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.grid()
ax1.set_xlabel("S (m)", fontsize=fontsize)
ax1.set_ylabel("Ti (eV)", fontsize=fontsize)
fig.tight_layout()
fig.show()
