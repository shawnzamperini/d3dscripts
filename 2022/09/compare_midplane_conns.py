# Script to compare midplane connection lengths from some MAFOT simulations.
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import numpy as np
from matplotlib.colors import LogNorm, Normalize


# Get the dictionary with the number of coordinates per flux tube.
ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/ncoords_167463_tor.pickle"
with open(ncoords_path, "rb") as f:
    ncoords = pickle.load(f)


# Load in the mafot data into a dataframe.
root = "/Users/zamperini/Documents/d3d_work/mafot_files/167463/"
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
def load_mafot(mafot_path):
    df = pd.read_csv(mafot_path, skiprows=52, names=columns, delimiter="\t")
    r = df["R (m)"].values
    z = df["Z (m)"].values
    l = df["Lconn (km)"].values * 1000

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

    return conns

polerr = load_mafot(root+"lam_pol.dat")
torerr = load_mafot(root+"lam_tor.dat")
polnoerr = load_mafot(root+"lam_polnoerr.dat")
tornoerr = load_mafot(root+"lam_tornoerr.dat")


psin_idx = 16  # About psin = 1.40.
def get_plot_params(conns):
    psins = list(conns.keys())
    psin = psins[psin_idx]
    l = conns[psin]["l"]
    d = conns[psin]["d"]
    z = conns[psin]["z"]

    # Set zeros to something outside the colorbar range so they don't show up.
    l[l==0] = -999

    return {"d":d, "l":l, "z":z, "psin":psin}


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8))

polerr_params = get_plot_params(polerr)
torerr_params = get_plot_params(torerr)
polnoerr_params = get_plot_params(polnoerr)
tornoerr_params = get_plot_params(tornoerr)
params = [polnoerr_params, tornoerr_params, polerr_params, torerr_params]

# For the colorbar.
minz = -1.5286
maxz = 1.5119
minl = 0.0
maxl = 13.623
norm = Normalize(minl, maxl)
levels = np.arange(0, maxl, 2)

for ax, p in zip([ax1, ax2, ax3, ax4], params):
    tri = ax.tricontourf(p["d"], p["z"], p["l"], cmap="inferno", norm=norm,
        levels=levels, extend="max")
    ax.set_facecolor("grey")
    ax.set_ylim(minz, maxz)
    psin_label = r"$\mathdefault{\psi_n}$" + " = {:.3f}".format(p["psin"])
    ax.text(0.5, 0.90, psin_label, transform=ax.transAxes,
        bbox={"facecolor":"white", "edgecolor":"k"}, fontsize=16)

ax1.set_ylabel("No Error Fields", fontsize=16)
ax1.set_title("Current Wall", fontsize=16)
ax2.set_title("Toroidal Limiters", fontsize=16)
ax3.set_ylabel("With Error Fields", fontsize=16)

fig.supxlabel("Toroidal Angle", fontsize=16)
fig.supylabel("Z (m)", fontsize=16)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cbar = fig.colorbar(tri, cax=cbar_ax)
cbar.set_label("Connection Length (m)", fontsize=16)

#fig.tight_layout()
fig.show()
