import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-bkg-001.nc"
op = oedge_plots.OedgePlots(ncpath)

# Get LFS neutral density. Aaron said the HFS isn't good data.
llama_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/LLAMA_190422.npz"
llama = np.load(llama_path)
lrho = llama["rho"]
lpsin = np.sqrt(lrho)  # According to Florian L.
lneut = llama["nDens_LFS"].mean(axis=0)
#lneut_err = np.sqrt(np.square(llama["nDens_LFS_err"]).sum(axis=0))
lneut_err = np.nanmean(llama["nDens_LFS_err"]mean(axis=0)

# Get fake probe of neutral density at the LLAMA location.
neut_probe = op.fake_probe(1.836, 2.07, -0.778, -0.76, data="neut_dens", plot="psin", show_plot=False)
neut_psin = neut_probe["psin"]

fig, ax = plt.subplots(figsize=(5,4))
#ax.fill_between(lpsin, lneut-lneut_err/2, lneut+lneut_err/2)
ax.plot(lpsin, lneut, label="LLAMA")
ax.plot(neut_psin, neut_probe["neut_dens"], label="OSM")
ax.legend()
ax.grid(which="both", alpha=0.5)
ax.set_xlim([0.9, 1.15])
ax.set_ylim([1e13, 5e16])
ax.set_yscale("log")
ax.set_xlabel("Psin")
ax.set_ylabel("Neutral Density (m-3)")
fig.tight_layout()
fig.show()
