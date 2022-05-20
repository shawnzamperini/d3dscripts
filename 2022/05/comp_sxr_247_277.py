# Compare the SXR data for 167247 and 167277.
import numpy as np
import matplotlib.pyplot as plt


path247 = "/Users/zamperini/Documents/d3d_work/files/imp_analysis_167247.npz"
path277 = "/Users/zamperini/Documents/d3d_work/files/imp_analysis_167277.npz"
sxr247 = np.load(path247)
sxr277 = np.load(path277)

rho247 = sxr247["rho"]
cw247 = sxr247["cw"]
cw247_avg = np.nanmean(cw247, axis=0)
cw247_std = np.nanstd(cw247, axis=0)
rho277 = sxr277["rho"]
cw277 = sxr277["cw"]
cw277_avg = np.nanmean(cw277, axis=0)
cw277_std = np.nanstd(cw277, axis=0)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(rho247, cw247_avg*100, color="k", label="167247")
ax.fill_between(rho247, (cw247_avg-cw247_std)*100, (cw247_avg+cw247_std)*100, color="k", alpha=0.3)
ax.plot(rho277, cw277_avg*100, color="r", label="167277")
ax.fill_between(rho277, (cw277_avg-cw277_std)*100, (cw277_avg+cw277_std)*100, color="r", alpha=0.3)
ax.legend(fontsize=14)
#ax.set_xlim([0.7, 1])
#ax.set_ylim([0, 0.5])
ax.set_xlabel(r"$\mathdefault{\rho}$", fontsize=14)
ax.set_ylabel("W Concentration (%)", fontsize=14)
fig.tight_layout()
fig.show()
