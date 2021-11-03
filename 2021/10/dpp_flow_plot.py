# Script to plot the flow for our given shot.
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np


# Path to RCP data.
rcp_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
rcp1 = pd.read_excel(rcp_path, sheet_name="MP187111_1")
rcp2 = pd.read_excel(rcp_path, sheet_name="MP187111_2")

rmrs1 = rcp1["Dist. from sep."].values * 100
rmrs2 = rcp2["Dist. from sep."].values * 100
mach1 = rcp1["Mach"].values
mach2 = rcp2["Mach"].values

# Let's just do an interpolation and average them.
f1 = interp1d(rmrs1, mach1)
f2 = interp1d(rmrs2, mach2)
#com_rmrs = np.linspace(max(rmrs1.min(), rmrs2.min()), min(rmrs1.max(), rmrs2.max()), 50)
com_rmrs = np.linspace(max(rmrs1.min(), rmrs2.min()), 9.0, 50)
avg_mach = (f1(com_rmrs) + f2(com_rmrs)) / 2.0

fig, ax = plt.subplots(figsize=(5,4))
#ax.plot(rmrs1, mach1)
#ax.plot(rmrs2, mach2)
ax.plot(com_rmrs, avg_mach, lw=3, color="tab:red")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim([-0.2, 0.2])
ax.set_xlabel(r"$\mathdefault{R-R_{sep}\ (m)}$", fontsize=14)
ax.set_ylabel("Mach Number", fontsize=14)
ax.tick_params(which="both", labelsize=12)
ax.axhline(0.0, color="k", lw=3)
ax.grid()
fig.tight_layout()
fig.show()
