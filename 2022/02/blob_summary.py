import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 5))

for shot in [167193, 167195, 184527, 187111]:
    x = dfs[shot]["R-Rsep(cm)"]
    vr = dfs[shot]["Vr(m/s)"]
    fblob = dfs[shot]["Npeaks"] / 5e-3  # 5ms intervals
    drad = dfs[shot]["D_rad(cm)"]
    epol = dfs[shot]["Epol(V/m)"]
    ax1.plot(x, vr, label=shot)
    ax2.plot(x, fblob, label=shot)
    ax3.plot(x, drad, label=shot)
    ax4.plot(x, epol, label=shot)

    print("{}: {:.2f} Hz".format(shot, fblob.mean()))

ax1.legend()
ax1.set_xlabel("R-Rsep (cm)")
ax1.set_ylabel("Blob velocity (m/s)")
ax2.set_xlabel("R-Rsep (cm)")
ax2.set_ylabel("Blob frequency (s-1)")
ax3.set_xlabel("R-Rsep (cm)")
ax3.set_ylabel("Blob width (cm)")
ax4.set_xlabel("R-Rsep (cm)")
ax4.set_ylabel("Epol (V/m)")

fig.tight_layout()
fig.show()
