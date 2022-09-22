# Script to compare DIVIMP flows at MiMES with the Mach probe measurements.
import oedge_plots
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Load measurements. Only care about 194 and 195 as those were the shapes we used.
# M > 0 = towards inner target?
xlpath = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_167192-195.xlsx"
plunges = ["MP167194_1", "MP167194_2", "MP167195_1", "MP167195_2"]
dfs = []
for plunge in plunges:
    dfs.append(pd.read_excel(xlpath, sheet_name=plunge))

# Load DIVIMP flows.
#ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div14a.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-bg-shifted-ring-entry-09.nc"
op = oedge_plots.OedgePlots(ncpath)
mdict = op.fake_probe(r_start=2.22, r_end=2.37, z_start=-0.188, z_end=-0.188, data="Mach", skip_t13=True, plot="R")
r_t13 = mdict["r"]

# For DIVIMP, M > 0 = towards outer target, so multiply by -1 to agree with RCP.
m_t13 = np.array(mdict["Mach"]) * -1

# Can't hurt to throw the DIVIMP connection length on top.
#lotfdict = op.fake_probe(2.20, 2.38, -0.185, -0.185, "L OTF")
#litfdict = op.fake_probe(2.20, 2.38, -0.185, -0.185, "L ITF")
#lotf = np.array(lotfdict["L OTF"])
#litf = np.array(litfdict["L ITF"])

fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.axhline(0, color="k")
for i in range(len(dfs)):
    df = dfs[i]
    r = df["R(cm)"] / 100
    m = df["Machn"]
    ax1.plot(r, m, color="k", alpha=0.4)
ax1.plot(r_t13, m_t13, color="r")
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Mach")
ax1.set_ylim(-1, 1)
#ax11 = ax1.twinx()
#ax11.plot(r_t13, lotf+litf, color="k", linestyle="--")
#ax11.set_yscale("log")
#ax11.set_ylabel("Connection Length (m)")
fig.tight_layout()
fig.show()
