# This script will compare the calculated radial profiles of the flow in our
# narrow grid to those measured form the RCP plunges.
import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge/")
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mach = True

ncpath = "/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-001.nc"
op = oedge_plots.OedgePlots(ncpath)

xl_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/" + \
  "rcp_master_detailed.xlsx"
master_xl = pd.ExcelFile(xl_path)


# Get the RCP data.
shot = 184267
mp_dfs = []
mp_dfs.append(master_xl.parse("MP{}_1".format(shot)))
mp_dfs.append(master_xl.parse("MP{}_2".format(shot)))
xp_dfs = []
xp_dfs.append(master_xl.parse("XP{}_1".format(shot)))
xp_dfs.append(master_xl.parse("XP{}_2".format(shot)))

if mach:
    rcp_data = "Mach"
    op_data = "Mach"
else:
    rcp_data = "Vflow (m/s)"
    op_data = "Velocity"

# Get the OSM "plunges".
op_mp_x, op_mp_y = op.fake_probe(2.25, 2.28, -0.188, -0.188, data=op_data,
  show_plot=False, plot="psin")
op_xp_x, op_xp_y = op.fake_probe(1.493, 1.493, -1.01, -1.08, data=op_data,
  show_plot=False, plot="psin")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))

ax1.plot(mp_dfs[0]["Psin"], mp_dfs[0][rcp_data], color="r")
ax1.plot(mp_dfs[1]["Psin"], mp_dfs[1][rcp_data], color="r")
ax1.plot(op_mp_x, op_mp_y, color="g")

ax2.plot(xp_dfs[0]["Psin"], xp_dfs[0][rcp_data], color="r", label="RCP")
ax2.plot(xp_dfs[1]["Psin"], xp_dfs[1][rcp_data], color="r")
ax2.plot(op_xp_x, op_xp_y, color="g", label="OSM")

ax1.set_xlabel("Psin")
ax1.set_ylabel(op_data)
ax2.set_xlabel("Psin")
ax1.set_title("MRCP")
ax2.set_title("XRCP")
if mach:
    ax1.set_ylim([-1, 1])
    ax2.set_ylim([-1, 1])
ax1.grid()
ax2.grid()
ax2.legend()

fig.tight_layout()
fig.show()
