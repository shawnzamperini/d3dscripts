# This script is similar to 03/conn_psin_gif.py, except it's a bit simpler
# in that it just plots the connection lengths at a specified Z location.
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

plt.rcParams['font.family'] = 'Century Gothic'
# Inputs
#shot = 176971
shot = 167196
pol_lims = False
include_coils = True

# Set the correct paths for each case.
gfile_pickle_path = "/Users/zamperini/Documents/d3d_work/archive/167196/167196_3500.pickle"
pol_wall_path = "/Users/zamperini/Documents/d3d_work/lwt/167196/mafot_3d_wall.dat"
tor_wall_path = "/Users/zamperini/Documents/d3d_work/lwt/930116/mafot_wall_wide.dat"
if shot == 167196:
    pol_mafot_path_fci = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_167196_with_fci_coils_parametrized.dat"
    pol_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_with_pol_lims.dat"
    pol_ncoords_path = "/Users/zamperini/github/d3dscripts/2022/03/ncoords_167196_pol.pickle"
    tor_mafot_path_fci = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_167196_with_fci_coils_parametrized_tl.dat"
    tor_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/lam_with_tor_lims.dat"
    tor_ncoords_path = "/Users/zamperini/github/d3dscripts/2022/03/ncoords_167196_tor.pickle"
    vmin = 0
    vmax = 20
elif shot == 186914:
    pol_mafot_path_fci = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/lam_parametrized_pol.dat"
    pol_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/lam_parametrized_pol.dat"
    pol_ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/ncoords_186914_pol.pickle"
    tor_mafot_path_fci = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/lam_parametrized_tor.dat"
    tor_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/lam_parametrized_tor.dat"
    tor_ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/186914/ncoords_186914_tor.pickle"
    vmin = 0
    vmax = 12
elif shot == 176971:
    pol_mafot_path_fci = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/lam_paramterized_pol.dat"
    pol_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/lam_paramterized_pol.dat"
    pol_ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/ncoords_176971_pol.pickle"
    tor_mafot_path_fci = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/lam_paramterized_tor.dat"
    tor_mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/lam_paramterized_tor.dat"
    tor_ncoords_path = "/Users/zamperini/Documents/d3d_work/mafot_files/176971/ncoords_176971_tor.pickle"
    vmin = 0
    vmax = 12

def load_conns(mafot_path, ncoords_path, wall_path):

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

    return conns

print("pol_fci")
pol_fci = load_conns(pol_mafot_path_fci, pol_ncoords_path, pol_wall_path)
print("pol")
pol = load_conns(pol_mafot_path, pol_ncoords_path, pol_wall_path)
print("tor_fci")
tor_fci = load_conns(tor_mafot_path_fci, tor_ncoords_path, tor_wall_path)
print("tor")
tor = load_conns(tor_mafot_path, tor_ncoords_path, tor_wall_path)



#print("plot_z = {:.4f} m".format(plot_z))
pol_ls = []
pol_fci_ls = []
tor_ls = []
tor_fci_ls = []
plot_psins_list = []
for ls, conns in zip([pol_ls, pol_fci_ls, tor_ls, tor_fci_ls], [pol, pol_fci, tor, tor_fci]):
    psins = np.array(list(conns.keys()))
    plot_psins = psins
    plot_psins_list.append(plot_psins)
    for psin in plot_psins:
        zs = np.unique(conns[psin]["z"])
        idx = np.argmin(np.abs(zs))
        mask = np.array(conns[psin]["z"]) == zs[idx]
        if mask.sum() == 0:
            print("Error! mask.sum() == 0")
        degs = np.array(conns[psin]["d"])[mask]
        rs = np.array(conns[psin]["r"])[mask]
        ls.append([degs, np.array(conns[psin]["l"])[mask], rs])

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))

for i in range(0, len(plot_psins)):
    d = pol_ls[i][0]
    l = pol_ls[i][1]
    ax1.plot(d, l, label="{:.2}".format(plot_psins[i]))
    d = pol_fci_ls[i][0]
    l = pol_fci_ls[i][1]
    ax2.plot(d, l)
    d = tor_ls[i][0]
    l = tor_ls[i][1]
    ax3.plot(d, l)
    d = tor_fci_ls[i][0]
    l = tor_fci_ls[i][1]
    ax4.plot(d, l)

ax1.set_title("Current")
ax2.set_title("Current w/ FCI Coils")
ax3.set_title("Toroidal Limiters")
ax4.set_title("Toroidal Limiters w/ FCI Coils")
ax1.legend()
fig.tight_layout()
fig.show()


# Create a grid for pcolormesh with degree on the X axis, psin on Y axis and
# the connection length as the Z data.
pol_data = []
pol_fci_data = []
tor_data = []
tor_fci_data = []
for data, conns, plot_psins in zip([pol_data, pol_fci_data, tor_data, tor_fci_data], [pol_ls, pol_fci_ls, tor_ls, tor_fci_ls], plot_psins_list):
    dim1 = len(pol_ls[0][0]) + 1  # Number of degrees
    dim2 = len(plot_psins)  # Number of psins

    ls2d = np.zeros((dim1, dim2))
    dg2d = np.zeros(ls2d.shape)
    ps2d = np.zeros(ls2d.shape)
    rs2d = np.zeros(ls2d.shape)
    for d1 in range(0, dim1):
        for d2 in range(0, dim2):
            if d1 < dim1-1:
                dg2d[d1][d2] = conns[d2][0][d1]
                ps2d[d1][d2] = plot_psins[d2]
                ls2d[d1][d2] = conns[d2][1][d1]
                rs2d[d1][d2] = conns[d2][2][d1]
            else:
                # Assign last entry the first to close the loop.
                dg2d[d1][d2] = dg2d[d1-1][d2] + 2
                ps2d[d1][d2] = ps2d[d1-1][d2]
                ls2d[d1][d2] = ls2d[d1-1][d2]
                rs2d[d1][d2] = rs2d[d1-1][d2]

    data.append(dg2d)
    data.append(ps2d)
    data.append(ls2d)
    data.append(rs2d)

cmap = "inferno"
fontsize = 16
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 6), sharey=True, sharex=True)
ax = fig.add_subplot(111, frameon=False)

vmin = np.min([pol_data[2], pol_fci_data[2], tor_data[2], tor_fci_data[2]])
vmax = np.max([pol_data[2], pol_fci_data[2], tor_data[2], tor_fci_data[2]])
#ax1.pcolormesh(pol_data[0], pol_data[1], pol_data[2], shading="auto")
#ax2.pcolormesh(pol_fci_data[0], pol_fci_data[1], pol_fci_data[2], shading="auto")
#ax3.pcolormesh(tor_data[0], tor_data[1], tor_data[2], shading="auto")
#ax4.pcolormesh(tor_fci_data[0], tor_fci_data[1], tor_fci_data[2], shading="auto")
cont1 = ax1.contourf(pol_data[0], pol_data[1], pol_data[2], vmin=vmin, vmax=vmax, cmap=cmap)
ax2.contourf(pol_fci_data[0], pol_fci_data[1], pol_fci_data[2], vmin=vmin, vmax=vmax, cmap=cmap)
ax3.contourf(tor_data[0], tor_data[1], tor_data[2], vmin=vmin, vmax=vmax, cmap=cmap)
ax4.contourf(tor_fci_data[0], tor_fci_data[1], tor_fci_data[2], vmin=vmin, vmax=vmax, cmap=cmap)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(cont1, cax=cbar_ax)
cbar.set_label("Connection Length (m)", fontsize=fontsize)

ax.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
ax.set_xlabel("Toroidal Angle", fontsize=fontsize)
ax.set_ylabel("Psin\n", fontsize=fontsize)
ax1.set_ylim([None, 1.265])
#plt.gca().invert_yaxis()
ax1.set_title("Current", fontsize=fontsize)
ax2.set_title("Current w/ FCI Coils", fontsize=fontsize)
ax3.set_title("Toroidal Limiters", fontsize=fontsize)
ax4.set_title("Toroidal Limiters w/ FCI Coils", fontsize=fontsize)
fig.suptitle("Connection Lengths at OMP", fontsize=fontsize)
#fig.tight_layout()
fig.show()


# A circular contour plot????
cmap = "winter"   # winter,
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), sharex=True, sharey=True, subplot_kw=dict(projection='polar'))
#fig, ax1 = plt.subplots(1, 1, figsize=(6, 6), subplot_kw=dict(projection='polar'))
ax = fig.add_subplot(111, frameon=False)
ax.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

cont1 = ax1.contourf(np.radians(pol_fci_data[0]), pol_fci_data[3], pol_fci_data[2], vmin=vmin, vmax=vmax, cmap=cmap)
ax2.contourf(np.radians(tor_fci_data[0]), tor_fci_data[3], tor_fci_data[2], vmin=vmin, vmax=vmax, cmap=cmap)

fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.05])
cbar = fig.colorbar(cont1, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Connection Length (m)", fontsize=fontsize)

yticks = [2.31, 2.34, 2.37]
ax1.set_yticks(yticks)
ax1.tick_params(axis="both", labelsize=12)
ax2.tick_params(axis="both", labelsize=12)
ax1.grid(False, axis="x")
ax1.grid(axis="y", alpha=0.5)
ax2.grid(False, axis="x")
ax2.grid(axis="y", alpha=0.5)
ax1.set_ylim([2.25, 2.375])
ax1.set_title("Current Limiters", fontsize=fontsize)
ax2.set_title("Toroidal Limiters", fontsize=fontsize)
#fig.tight_layout()
fig.show()
#plt.savefig("/Users/zamperini/My Drive/Research/Documents/2022/04/polar_limiter_comparison_fci.png", transparent=True)

# Same plot just individual.
fig, ax1 = plt.subplots(1, 1, figsize=(6, 6), subplot_kw=dict(projection='polar'))
ax = fig.add_subplot(111, frameon=False)
ax.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
cont1 = ax1.contourf(np.radians(pol_fci_data[0]), pol_fci_data[3], pol_fci_data[2], vmin=vmin, vmax=vmax, cmap=cmap)
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.05])
cbar = fig.colorbar(cont1, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Connection Length (m)", fontsize=fontsize)
yticks = [2.31, 2.34, 2.37]
ax1.set_yticks(yticks)
ax1.set_yticklabels(["R=2.31", "2.34", "2.37"])
ax1.tick_params(axis="both", labelsize=12)
ax1.grid(False, axis="x")
ax1.grid(axis="y", alpha=0.5)
ax1.set_ylim([2.25, 2.375])
ax1.set_title("Current Limiters", fontsize=fontsize)
#fig.tight_layout()
fig.show()
