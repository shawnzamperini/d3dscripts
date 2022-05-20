import numpy as np
import matplotlib.pyplot as plt


path = "/Users/zamperini/Documents/d3d_work/divimp_files/imp_analysis_167196.npz"
imps = np.load(path)
rho = imps["rho"]
cw = imps["cw"]

cw_avg = np.nanmean(cw, axis=0)
cw_std = np.nanstd(cw, axis=0)

fig, ax1 = plt.subplots()

ax1.plot(rho, cw_avg, color="r")
ax1.fill_between(rho, cw_avg-cw_std, cw_avg+cw_std, color="r", alpha=0.3)

ax1.set_xlabel("Rho")
ax1.set_ylabel("cW")

ax1.set_title("167196 2500-5000ms")

fig.tight_layout()
fig.show()
