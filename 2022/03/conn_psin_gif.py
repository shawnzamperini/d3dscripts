# This script works in tandem with paramtrize_psin.py, with MAFOT as the
# middleman. That script gives hthe points at which to run MAFOT for each
# flux surface, and this one takes the results and makes a gif out of them as
# one progress outwards along the flux surfaces.
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
import imageio
from tqdm import tqdm
import numpy as np
from matplotlib.colors import LogNorm, Normalize
import matplotlib.tri as tri
from LimWallToolkit import LimWallToolkit
import sys


# Inputs
shot = 174783
pol_lims = False
include_coils = True

# Set the correct paths for each case.
if pol_lims:
    wall_path = "/Users/zamperini/Documents/d3d_work/lwt/167196/mafot_3d_wall.dat"
    if shot == 167196:
        if include_coils:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_167196_with_fci_coils_parametrized.dat"
            gif_path = "psin_conns_167196_pol_coils.gif"
        else:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_with_pol_lims.dat"
            gif_path = "psin_conns_167196_pol.gif"
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/archive/167196/167196_3500.pickle"
        ncoords_path = "/Users/zamperini/github/d3dscripts/2022/03/ncoords_167196_pol.pickle"

    elif shot == 167463:
        if include_coils:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/lam_pol.dat"
            gif_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/psin_conns_167463_pol.gif"
        else:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/lam_polnoerr.dat"
            gif_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/psin_conns_167463_pol_noerr.gif"
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/167463_3000.pickle"
        ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/ncoords_167463_pol.pickle"

    elif shot == 174783:
        if include_coils:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/lam_svr_pol.dat"
            gif_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/psin_conns_178783_svr_pol.gif"
        else:
            print("Not run without coils yet.")
            sys.exit()
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/174783_0.pickle"
        ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/ncoords_174783_svr.pickle"

else:
    wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
    if shot == 167196:
        if include_coils:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_167196_with_fci_coils_parametrized_tl.dat"
            gif_path = "psin_conns_167196_tor_coils.gif"
        else:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_with_tor_lims.dat"
            gif_path = "psin_conns_167196_tor.gif"
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/archive/167196/167196_3500.pickle"
        ncoords_path = "/Users/zamperini/github/d3dscripts/2022/03/ncoords_167196_tor.pickle"

    elif shot == 167463:
        if include_coils:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/lam_tor.dat"
            gif_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/psin_conns_167463_tor.gif"
        else:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/lam_tornoerr.dat"
            gif_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/psin_conns_167463_tor_noerr.gif"
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/167463_3000.pickle"
        ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/ncoords_167463_tor.pickle"

    elif shot == 174783:
        if include_coils:
            mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/lam_svr_tor.dat"
            gif_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/psin_conns_178783_svr_tor.gif"
        else:
            print("Not run without coils yet.")
            sys.exit()
        gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/174783_0.pickle"
        ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/174783/ncoords_174783_svr.pickle"



# Load in the mafot data into a dataframe.
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
df = pd.read_csv(mafot_path, skiprows=52, names=columns, delimiter="\t")
r = df["R (m)"].values
z = df["Z (m)"].values
l = df["Lconn (km)"].values * 1000

# Get the dictionary with the number of coordinates per flux tube.
with open(ncoords_path, "rb") as f:
    ncoords = pickle.load(f)

# Load in the wall coordinates and gfile.
lwt = LimWallToolkit()
wall = lwt.read_3d_wall(wall_path)

with open(gfile_pickle_path, "rb") as f:
    gfile = pickle.load(f)
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin_rz = gfile["PSIRZ_NORM"]
R_sep = gfile["RBBBS"]
Z_sep = gfile["ZBBBS"]

# Limit to just the right of the OMP.
raxis = gfile["RMAXIS"]
outboard = gfile["R"] > raxis
gR_out = gR[:, outboard]
gZ_out = gZ[:, outboard]
psin_rz_out = psin_rz[:, outboard]

# The coordinates are printed one psin at a time, [0-360] before moving onto
# the next psin.
conns = {}
n = 0
minz = 999; maxz = -999
minl = 999; maxl = -999
for psin, ncoord in ncoords.items():
    tmp_l = []; tmp_r = []; tmp_z = []; tmp_d = []
    for deg in range(0, 360, 2):
        for i in range(0, ncoord):
            tmp_l.append(l[n])
            tmp_r.append(r[n])
            tmp_z.append(z[n])
            tmp_d.append(deg)

            # Keep track of minimum and maximum values for plots.
            if z[n] < minz:
                minz = z[n]
            if z[n] > maxz:
                maxz = z[n]
            if l[n] < minl:
                minl = l[n]
            if l[n] > maxl:
                maxl = l[n]

            n += 1
    conns[psin] = {"l":np.array(tmp_l), "r":np.array(tmp_r),
        "z":np.array(tmp_z), "d":np.array(tmp_d)}

# These should match, otherwise we missed a point.
if len(df) != n:
    print("Warning! Length of Dataframe and read in coordinates do no match!")
    print("len(df) = {}".format(len(df)))
    print("n = {}".format(n))

# For the colorbar.
norm = Normalize(minl, maxl)
levels = np.arange(0, maxl, 2)

# Plotting loop.
fnames = []
psins = list(conns.keys())
for i in tqdm(range(0, len(psins))):

    l = conns[psins[i]]["l"]
    d = conns[psins[i]]["d"]
    z = conns[psins[i]]["z"]

    # Set zeros to something outside the colorbar range so they don't show up.
    l[l==0] = -999

    fig, (ax2, ax) = plt.subplots(1, 2, figsize=(9, 4))
    tri = ax.tricontourf(d, z, l, cmap="inferno", norm=norm, levels=levels,
        extend="max")
    ax.set_facecolor("grey")
    ax.set_ylim(minz, maxz)
    psin_label = r"$\mathdefault{\psi_n}$" + " = {:.3f}".format(psins[i])
    ax.text(0.5, 0.90, psin_label, transform=ax.transAxes,
        bbox={"facecolor":"white", "edgecolor":"k"}, fontsize=16)
    ax.set_xlabel("Toroidal Angle", fontsize=16)
    ax.set_ylabel("Z (m)", fontsize=16)
    cbar = fig.colorbar(tri, ax=ax, ticks=levels)
    cbar.set_label("Connection Length (m)", fontsize=16)

    # Include a cross section that shows the flux surface as we are progressing
    flux_bounds = [psins[0], psins[-1]]
    ax2.contour(gR, gZ, psin_rz, levels=flux_bounds, colors="k", linewidths=1)
    ax2.contour(gR, gZ, psin_rz, levels=[1.0], linewidths=3, colors="k")
    ax2.contour(gR_out, gZ_out, psin_rz_out, levels=[psins[i]], colors="tab:red")
    ax2.plot(wall[242][0], wall[242][1], lw=2, color="k")
    ax2.set_aspect("equal")
    ax2.axis("off")

    fig.tight_layout()
    fname = "mafot_plots/plot_{}.png".format(i)
    fnames.append(fname)
    fig.savefig(fname)
    plt.close(fig)

# Now build the gif and clean up by removing the plots.
print("Saving as gif...")
with imageio.get_writer(gif_path, mode="I") as writer:
    for fname in fnames:
        image = imageio.imread(fname)
        writer.append_data(image)
for fname in fnames:
    os.remove(fname)
