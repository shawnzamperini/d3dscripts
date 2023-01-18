import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


path = "/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/blob_analyzed.xlsx"
xl = pd.ExcelFile(path)
sheets = xl.sheet_names

dfs = {}
for sheet in sheets:
    dfs[sheet] = xl.parse(sheet)

    # Data had arced.
    if sheet == "167195_1":
        dfs[sheet] = dfs[sheet].iloc[:4]

# Join like shots together.
dfs[167193] = dfs["167193_1"].append(dfs["167193_2"]).sort_values("R-Rsep(cm)")
dfs[167195] = dfs["167195_1"].append(dfs["167195_2"]).sort_values("R-Rsep(cm)")
dfs[184527] = dfs["184527_1"].append(dfs["184527_2"]).sort_values("R-Rsep(cm)")
dfs[187111] = dfs["187111_1"].append(dfs["187111_2"]).sort_values("R-Rsep(cm)")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5))

ax1.plot(dfs[167195]["R-Rsep(cm)"], dfs[167195]["Vr(m/s)"], lw=3, color="tab:red")
#ax1.plot([2, 5], [863, 863], linestyle="--", color="tab:red", lw=3)
ax1.plot([7, 10], [500, 500], linestyle="--", color="k", lw=3)
#ax1.set_xlabel(r"R-$\mathdefault{R_{sep}}$ (cm)", fontsize=14)
ax1.set_ylabel("Radial Blob\nVelocity (m/s)", fontsize=14)
ax1.grid()
ax1.set_ylim([0, 1200])

timestep = 0.005
fblob = dfs[167195]["Npeaks"] / timestep
ax2.plot(dfs[167195]["R-Rsep(cm)"], fblob, lw=3, color="tab:red")
ax2.set_ylabel("Blob Frequency", fontsize=14)
ax2.set_xlabel(r"R-$\mathdefault{R_{sep}}$ (cm)", fontsize=14)
ax2.grid()

fig.tight_layout()
fig.show()