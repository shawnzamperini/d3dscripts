# This script makes a gif of the impurity density as one moves outward.
import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge")
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import imageio
import os
from scipy.interpolate import interp1d


extra_frames = 10
num_s = 150
norm = False

ncpath = "/Users/zamperini/Documents/d3d_work/184527/d3d-184527-inj-007.nc"
op = oedge_plots.OedgePlots(ncpath)

# First load all the along ring data.
vals = []; psins = []
for ir in range(13, 38):
    #if ir in [23, 24, 31]:
    #    continue
    s, imp = op.along_ring(ir, "DDLIMS", charge="all", plot_it=False)

    # Want to normalize S.
    snorm = s / s.max()

    # Also want to inteprolate onto a common amount of X values for later.
    s_com = np.linspace(0, 1, 100)
    f_imp = interp1d(snorm, imp, fill_value="extrapolate")
    imp_com = f_imp(s_com)

    if norm:
        imp_com = imp_com / imp_com.max()

    vals.append([s_com, imp_com, s, imp])

    # Get the psin at this ring.
    _, psin = op.along_ring(ir, "PSIFL", plot_it=False)
    psin = np.array(psin)
    psins.append(psin[psin != 0].mean())

fnames = []
for ir in range(0, len(vals)-1):
    print("Plot {}...".format(ir))
    if ir in [22-13, 23-13, 24-13, 31-13]:
        continue

    snorm_diff = vals[ir+1][0] - vals[ir][0]
    imp_diff = vals[ir+1][1] - vals[ir][1]
    psin_diff = psins[ir+1] - psins[ir]
    for j in range(0, extra_frames+1):
        interp_snorm = vals[ir][0] + (snorm_diff / extra_frames) * j
        interp_imp = vals[ir][1] + (imp_diff / extra_frames) * j
        interp_psin = psins[ir] + (psin_diff / extra_frames) * j

        fig, ax = plt.subplots()

        label = r"$\phi_n$" + " = {:.3f}".format(interp_psin)
        #label = "ir = {}".format(ir)
        ax.plot(interp_snorm, interp_imp, color="tab:red", lw=2, label=label)
        ax.legend(loc="upper right")
        ax.set_xlabel("Distance from outer target (normalized)", fontsize=14)
        ax.set_ylabel("C13 Density (m-3)", fontsize=14)

        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1e16])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid()

        fig.tight_layout()

        # Save figure and close.
        fname = "sep_out_plots/{:}_{:}.png".format(ir,j)
        fnames.append(fname)
        fig.savefig(fname)
        plt.close(fig)

# Now build the gif and clean up by removing the plots.
print("Saving as gif...")
with imageio.get_writer('sep_outwards.gif', mode="I") as writer:
    for fname in fnames:
        image = imageio.imread(fname)
        writer.append_data(image)
for fname in fnames:
    os.remove(fname)
