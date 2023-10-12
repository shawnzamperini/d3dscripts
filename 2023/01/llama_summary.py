# Plot together the LLAMA data from 190484-6.
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt


# Note the underscore at the end. Aaron added them to indicate they have been re-run after fixing a bug.
llama4 = np.load("/Users/zamperini/Documents/d3d_work/files/LLAMA_190484_.npz")
llama5 = np.load("/Users/zamperini/Documents/d3d_work/files/LLAMA_190485_.npz")
llama6 = np.load("/Users/zamperini/Documents/d3d_work/files/LLAMA_190486_.npz")

# 190485 is a different file format than the other two.
# lpsin5 = np.tile(np.square(llama5["rho"]), 98)
# lneut5 = llama5["nDens_LFS"].flatten()
# lneut_err5 = llama5["nDens_LFS_err"].flatten()
# lion5 = llama5["ion_LFS"].flatten()
# lion_err5 = llama5["ion_LFS_err"].flatten()

# Now do 484 and 486.
lpsin4 = llama4["psi_n"]
lpsin5 = llama5["psi_n"]
lpsin6 = llama6["psi_n"]
lneut4 = llama4["nDens_LFS"]
lneut5 = llama5["nDens_LFS"]
lneut6 = llama6["nDens_LFS"]
lneut_err4 = llama4["nDens_LFS_err"]
lneut_err5 = llama5["nDens_LFS_err"]
lneut_err6 = llama6["nDens_LFS_err"]
lion4 = llama4["ion_LFS"]
lion5 = llama5["ion_LFS"]
lion6 = llama6["ion_LFS"]
lion_err4 = llama4["ion_LFS_err"]
lion_err5 = llama5["ion_LFS_err"]
lion_err6 = llama6["ion_LFS_err"]

# Some basic filtering.
filt = medfilt
window = 25
lion4_filt = filt(lion4, window)
lion5_filt = filt(lion5, window)
lion6_filt = filt(lion6, window)
lneut4_filt = filt(lneut4, window)
lneut5_filt = filt(lneut5, window)
lneut6_filt = filt(lneut6, window)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharex=True)

ax1.axvline(1.0, color="k", zorder=5)
ax1.scatter(lpsin4, lion4, label="190484", s=15, edgecolors="k", color="tab:red", zorder=15)
ax1.scatter(lpsin5, lion5, label="190485", s=15, edgecolors="k", color="tab:purple", zorder=15)
ax1.scatter(lpsin6, lion6, label="190486", s=15, edgecolors="k", color="tab:cyan", zorder=15)
ax1.plot(lpsin4, lion4_filt, color="k", lw=3, zorder=25)
ax1.plot(lpsin5, lion5_filt, color="k", lw=3, zorder=25)
ax1.plot(lpsin6, lion6_filt, color="k", lw=3, zorder=25)
ax1.plot(lpsin4, lion4_filt, color="tab:red", lw=2, zorder=25)
ax1.plot(lpsin5, lion5_filt, color="tab:purple", lw=2, zorder=25)
ax1.plot(lpsin6, lion6_filt, color="tab:cyan", lw=2, zorder=25)
#ax1.legend(fontsize=14)
ax1.set_yscale("log")
ax1.grid(which="both", alpha=0.3)
ax1.set_xlim([0.95, 1.10])
ax1.set_ylim([1e19, 1e21])
ax1.set_xlabel("Psin", fontsize=14)
ax1.set_ylabel("Ionization Rate (m-3 s-1)", fontsize=14)

ax2.axvline(1.0, color="k", zorder=5)
ax2.scatter(lpsin4, lneut4, label="190484", s=15, edgecolors="k", color="tab:red", zorder=15)
ax2.scatter(lpsin5, lneut5, label="190485", s=15, edgecolors="k", color="tab:purple", zorder=15)
ax2.scatter(lpsin6, lneut6, label="190486", s=15, edgecolors="k", color="tab:cyan", zorder=15)
ax2.plot(lpsin4, lneut4_filt, color="k", lw=3, zorder=25)
ax2.plot(lpsin5, lneut5_filt, color="k", lw=3, zorder=25)
ax2.plot(lpsin6, lneut6_filt, color="k", lw=3, zorder=25)
ax2.plot(lpsin4, lneut4_filt, color="tab:red", lw=2, zorder=25)
ax2.plot(lpsin5, lneut5_filt, color="tab:purple", lw=2, zorder=25)
ax2.plot(lpsin6, lneut6_filt, color="tab:cyan", lw=2, zorder=25)
ax2.legend(fontsize=12)
ax2.set_yscale("log")
ax2.grid(which="both", alpha=0.3)
ax2.set_ylim([1e13, 1e18])
ax2.set_xlabel("Psin", fontsize=14)
ax2.set_ylabel("Neutral Density (m-3)", fontsize=14)

fig.tight_layout()
fig.show()