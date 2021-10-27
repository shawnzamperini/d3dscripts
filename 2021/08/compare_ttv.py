# This script is meant to create comparisons of the C2 and C3 TTV data with
# output from a DIVIMP run. The TTV data file convention is par = C2 and
# perp = C3. f120_f135 is from before the methane puff (times vary) and
# f234_f246 is from after the puff (times vary). The difference between these
# two give you the emission that is from the methane, in theory.
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from scipy.interpolate import interp2d
import matplotlib.colors as colors
import matplotlib as mpl

import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge/")
import oedge_plots
sys.path.append("/Users/zamperini/github/utk-fusion/tools/")
from load_gfile_mds import load_gfile_mds


# Some constants.
shot = 184527
efit = "EFIT_CAKE02"
#efit = "EFIT02"
ncpath = "/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-001.nc"
sav_base = "/Users/zamperini/Documents/d3d_work/methane_ttv_data/emission" + \
  "_structure_out_cam0"

def load_avg_inverted(path, gfile):
    """
    Load in the sav file, and then get the average inverted profile for each
    time (axis 0).
    """
    sav = readsav(path)
    inverted = sav("emission_structure")["inverted"][0].mean(axis=0)
    r = sav("emission_structure")["radii"][0][0]
    z = sav("emission_structure")["elevation"][0][0]
    R, Z = np.meshgrid(r, z)

    # Load the gfile so that we can map the TTV (R, Z) to psin.
    f_psin = interp2d(gfile["R"], gfile["Z"], gfile["psiRZn"])
    psin = f_psin(r, z)

    return {"R":R, "Z":Z, "psin":psin, "inverted":inverted}

# Load in everything we need.
op = oedge_plots.OedgePlots(ncpath)
gfile = load_gfile_mds(shot=shot, time=3000, tree=efit, tunnel=False)
c2_before = load_avg_inverted(sav_base + "par_"  + str(shot) + "_f102_f114.sav", gfile)
c2_after  = load_avg_inverted(sav_base + "par_"  + str(shot) + "_f222_f234.sav", gfile)
c3_before = load_avg_inverted(sav_base + "perp_" + str(shot) + "_f102_f114.sav", gfile)
c3_after  = load_avg_inverted(sav_base + "perp_" + str(shot) + "_f222_f234.sav", gfile)

c2_diff_inv = c2_after["inverted"] - c2_before["inverted"]
c3_diff_inv = c3_after["inverted"] - c3_before["inverted"]

# Setup a symlog colorbar.
cmap = "coolwarm"
vmin = c2_diff_inv.min()
vmax = c2_diff_inv.max()
norm = mpl.colors.SymLogNorm(linthresh=0.01 * vmax, vmin=vmin, vmax=vmax, base=10)
scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

#ax.contourf(c2_before["R"], c2_before["Z"], c2_before["inverted"])

# Figure out where the colorbar labels will be. Here, log scale.
cbar_max = 0.3
levels = np.append(-np.geomspace(0.01, cbar_max, 5)[::-1],
  np.geomspace(0.01, cbar_max, 5))

# Complicated command to plot correctly.
cont1 = ax1.contourf(c2_before["R"], c2_before["Z"], c2_diff_inv,
  norm=colors.SymLogNorm(linthresh=0.005, vmin=-cbar_max, vmax=cbar_max, base=10),
  vmin=-cbar_max, vmax=cbar_max, cmap="coolwarm", levels=levels, extend="both")
cbar1 = fig.colorbar(cont1, ax=ax1, format="%.2f")
cbar1.ax.set_yticklabels(np.round(levels, 2))

cont2 = ax2.contourf(c3_before["R"], c3_before["Z"], c3_diff_inv,
  norm=colors.SymLogNorm(linthresh=0.01, vmin=-cbar_max, vmax=cbar_max, base=10),
  vmin=-cbar_max, vmax=cbar_max, cmap="coolwarm", levels=levels, extend="both")
cbar2 = fig.colorbar(cont2, ax=ax2, format="%.2f")
cbar2.ax.set_yticklabels(np.round(levels, 2))

# Equal aspect ratio.
ax1.axis("equal")
ax2.axis("equal")

# Plot wall and flux contours.
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
ax1.plot(gfile["wall"][:,0], gfile["wall"][:,1], color="k")
ax1.contour(gR, gZ, gfile["psiRZn"], levels=np.arange(0.95, 1.1, 0.01), colors="k", alpha=0.7)
ax2.plot(gfile["wall"][:,0], gfile["wall"][:,1], color="k")
ax2.contour(gR, gZ, gfile["psiRZn"], levels=np.arange(0.95, 1.1, 0.01), colors="k", alpha=0.7)

ax1.set_title("C2 Puff Difference")
ax2.set_title("C3 Puff Difference")

ax1.set_xlim([1.065, 1.9])
ax1.set_ylim([0.6, 1.2])
ax2.set_xlim([1.065, 1.9])
ax2.set_ylim([0.6, 1.2])
fig.tight_layout()
fig.show()
