import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import os
import imageio
from mafot_read_3d_wall import read_wall
from tqdm import tqdm


# Options for plotting.
show_points = False
max_tor = 10  # Inclusive, eventually leave at 359.
extra_frames = 1
plot_ratio = False
plot_diff  = True

wall_path = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"
root = "/Users/zamperini/Documents/d3d_work/184527/mafot/"
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin", "psimax",
  "psiav", "pitch angle", "yaw angle", "theta", "psi"]

# Read in wall as dictionary of points at each angle.
wall = read_wall(wall_path)

# Need to load all the data first before getting into plotting.
Ls = {}
first_run = True
for tor in range(0, max_tor+1):

    # Load data into DataFrame.
    fname = root + "lam_torang{}.dat".format(tor)
    df = pd.read_csv(fname, skiprows=52, names=columns, delimiter="\t")

    # Reshape into plottable arrays. The R, Z locations don't change among the
    # runs, so just load and shape them once.
    if first_run:
        r = df["R (m)"].unique()
        z = df["Z (m)"].unique()
        R, Z = np.meshgrid(r, z)
        first_run = False
    l = df["Lconn (km)"].values * 1000  # km to m
    L = l.reshape(len(r), len(z))
    Ls[tor] = L

# The changes our data into ratios of how much L changes wrt to tor = 0 degrees.
if plot_ratio:
    denom = np.copy(Ls[0])
    for tor in range(0, max_tor+1):
        Ls[tor] = Ls[tor] / denom
    vmin = 1
    vmax = 2
    levels = np.linspace(vmin, vmax, 15)
    norm = Normalize(vmin, vmax)
    cbar_label = r"L / $L_{0^{\circ}}$"
    cmap = "Reds"
elif plot_diff:
    denom = np.copy(Ls[0])
    for tor in range(0, max_tor+1):
        Ls[tor] = Ls[tor] - denom
    vmin = -10
    vmax = 10
    levels = np.linspace(vmin, vmax, 15)
    norm = Normalize(vmin, vmax)
    cbar_label = r"L - $\mathdefault{L_{0^{\circ}}}$"
    cmap = "coolwarm"
else:
    vmin = 0.1
    vmax = 100
    levels = np.geomspace(vmin, vmax, 15)
    norm = LogNorm(vmin=vmin, vmax=vmax)
    cbar_label = "L"
    cmap = "Reds"

# Plotting loop.
fnames = []
for tor in tqdm(range(0, max_tor)):
    #print("Plot {}...".format(tor))

    L_diff = Ls[tor + 1] - Ls[tor]
    for j in range(0, extra_frames+1):

        interp_L = Ls[tor] + (L_diff / extra_frames) * j

        # To just silence a warning.
        interp_L = np.ma.masked_where(interp_L <= 0, interp_L)

        fig, ax = plt.subplots(figsize=(4, 7))
        #cont = ax.contourf(R, Z, interp_L,
        #  norm=norm, extend="both",
        #  cmap=cmap, levels=levels)
        cont = ax.pcolormesh(R, Z, interp_L,
          norm=norm,
          cmap=cmap)
        cbar = fig.colorbar(cont, shrink=0.5, ticks=levels, format="%3.1f")
        cbar.set_label(cbar_label, fontsize=14)

        if show_points:
            ax.scatter(R.flatten(), Z.flatten(), s=1, color="k")

        # 360 - angle because...
        ax.plot(wall[360-tor][0], wall[360-tor][1], color="k", lw=2)

        label = r"$\phi$ = " + str(tor) + r"$^\circ$"
        ax.text(0.6, 0.95, label, transform=ax.transAxes, fontsize=14)

        ax.set_xlim(0.95, 2.55)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect("equal")
        ax.set_xlabel("R (m)")
        ax.set_ylabel("Z (m)")
        fig.tight_layout()

        fname = "mafot_plots/plot_{}-{}.png".format(tor, j)
        fnames.append(fname)
        fig.savefig(fname)
        plt.close(fig)

# Now build the gif and clean up by removing the plots.
print("Saving as gif...")
with imageio.get_writer('mafot_2d.gif', mode="I") as writer:
    for fname in fnames:
        image = imageio.imread(fname)
        writer.append_data(image)
for fname in fnames:
    os.remove(fname)
