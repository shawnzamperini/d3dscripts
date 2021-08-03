# This script will create a 2D cross section of the impurity density profile
# and also highlight the far-SOL area that guides 3DLIM input.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
import matplotlib as mpl

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Get case representative of each direction.
unf_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031d.nc"
unf = oedge_plots.OedgePlots(unf_path)
fav = oedge_plots.OedgePlots(fav_path)

# Better wall data.
wall_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/grid-files/d3d_wall_mrc_mod.xlsx"
wall_df = pd.read_excel(wall_path, usecols=["rnew", "znew"]) / 1000

# Grab the 2D impurity data and normalize it.
unf_imp_data = unf.read_data_2d("DDLIMS", charge="all")
fav_imp_data = fav.read_data_2d("DDLIMS", charge="all")
max_imp = np.max((unf_imp_data, fav_imp_data))
unf_imp_data = unf_imp_data / max_imp
fav_imp_data = fav_imp_data / max_imp

# Use the unfavorable data.
fig = unf.plot_contour_polygon("KTEBS", own_data=unf_imp_data, cmap="magma",
  wall_data=(wall_df["rnew"], wall_df["znew"]), vmin=0, vmax=0.8)

# Pull out the axes object so we can add shapes to it.
ax = fig.axes[0]
fig.axes[1].remove()
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="both", which="both", labelleft=False, labelbottom=False, left=False, bottom=False)

# Colors from magma colormap.
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))

# Add CP.
cp_r_tip = 2.25
cp_z = -0.188
cp_width = 0.03
rect_anchor = (cp_r_tip, cp_z - cp_width / 2)  # Z location of anchor just the bottom of CP
rect = patches.Rectangle(rect_anchor, 2.3755-cp_r_tip, cp_width, ec="k", color=colors[2])
ax.add_patch(rect)

# First create a new mesh of only the farsol ring. Copied from init of oedge_plots.
farsol_ring = 70
mesh = []
num_cells = 0
for ir in range(unf.nrs):
    for ik in range(unf.nks[ir]):
        index = unf.korpg[ir,ik] - 1
        if unf.area[ir,ik] != 0.0:

            # Only append if on the farsol ring.
            if ir == farsol_ring:
                vertices = list(zip(unf.rvertp[index][0:4], unf.zvertp[index][0:4]))
                mesh.append(vertices)
                num_cells = num_cells + 1

# Create array full of 1's to fill in on the ring.
blanks = np.full(np.array(mesh).shape[0], 1)

# Add the collection to our figure.
coll = mpl.collections.PolyCollection(mesh, array=blanks, cmap="winter",
                                      #cmap=scalar_map.cmap,
                                      #norm=scalar_map.norm,
                                      edgecolors='none')
ax.add_collection(coll)
