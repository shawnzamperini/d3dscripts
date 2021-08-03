import matplotlib.pyplot as plt
import numpy as np
import oedge_plots
import pandas as pd
import matplotlib.patches as patches


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

ncpath = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
op = oedge_plots.OedgePlots(ncpath)

# TS data for later.
ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/ts167247_final_v2.xlsx"
ts_df = pd.read_excel(ts_path)

# Pull out the wall data and fix where the antenna on the outboard side juts
# in a little.
wall_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/grid-files/d3d_wall_mrc_mod.xlsx"
wall_df = pd.read_excel(wall_path, usecols=["rnew", "znew"]) / 1000
keep_idx = np.where(np.logical_and(op.rvesm[0] != 0, op.zvesm[0] != 0))
rvesm = np.append(op.rvesm[0][keep_idx], op.rvesm[0][keep_idx][0])
zvesm = np.append(op.zvesm[0][keep_idx], op.zvesm[0][keep_idx][0])
rvesm[62:72] = rvesm[62]
zvesm[63] = zvesm[62]
zvesm[70] = zvesm[62]

# Create data of a bunch of zeros and make a 2D plot with a colormap that will
# show nothing except the separatrix.
blank_data = np.zeros(op.read_data_2d("KTEBS").shape)
fig = op.plot_contour_polygon("KTEBS", own_data=blank_data, cmap="bwr",
  wall_data=(wall_df["rnew"], wall_df["znew"]), show_mr=True)

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
fontsize = 16

# Add CP.
cp_r_tip = 2.25
cp_z = -0.188
cp_width = 0.03
rect_anchor = (cp_r_tip, cp_z - cp_width / 2)  # Z location of anchor just the bottom of CP
#rect = patches.Rectangle(rect_anchor, 2.377-cp_r_tip, cp_width, ec="k", color=colors[2])
rect = patches.Rectangle(rect_anchor, 2.3755-cp_r_tip, cp_width, ec="k", color=colors[2])
ax.add_patch(rect)

# Add CP annotation.
ax.annotate("Collector\nProbe", (cp_r_tip, rect_anchor[1] + cp_width / 2),
  xytext=(1.9, -0.4), arrowprops=dict(facecolor="black", arrowstyle="-"),
  fontsize=fontsize, ha="center")

# Metal tiles. Do it here so we can make them a bit thicker and different color.
tile_height = 0.02
facecolor   = colors[4]
edgecolor   = 'black'
floor_xy    = (1.32, -1.363 - tile_height)
floor_width = 1.37 - 1.32
shelf_xy    = (1.404, -1.250 - tile_height)
shelf_width = 1.454 - 1.404
floor_rect  = patches.Rectangle(floor_xy, width = floor_width,
                                   height=tile_height,
                                   facecolor=facecolor,
                                   edgecolor=edgecolor)
shelf_rect  = patches.Rectangle(shelf_xy, width = shelf_width,
                                   height=tile_height,
                                   facecolor=facecolor,
                                   edgecolor=edgecolor)
ax.add_patch(floor_rect)
ax.add_patch(shelf_rect)

# Add metal ring annotation.
ax.annotate("Metal Rings", (floor_xy[0] + floor_width / 2, floor_xy[1]), xytext=(1.6, -1.5),
  arrowprops=dict(facecolor="black", arrowstyle="-"),
  fontsize=fontsize)
ax.annotate("           ", (shelf_xy[0] + shelf_width / 2, shelf_xy[1]), xytext=(1.6, -1.5),
  arrowprops=dict(facecolor="black", arrowstyle="-"),
  fontsize=fontsize)

# Coordinates for Thomson scattering. Get them from the TS file.
ts_core = ts_df[ts_df["System"] == "core"]
ts_zs = ts_core["Z (m)"].unique()
ts_rs = np.full(len(ts_zs), ts_core["R (m)"].unique()[0])
ax.scatter(ts_rs, ts_zs, s=3, color=colors[0])
ax.annotate("Thomson\nScattering", (ts_rs[0], ts_zs[15]), xytext=(1.4, -0.4),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

ts_div = ts_df[ts_df["System"] == "divertor"]
ts_zs = ts_div["Z (m)"].unique()
ts_rs = np.full(len(ts_zs), ts_div["R (m)"].unique()[0])
ax.scatter(ts_rs, ts_zs, s=3, color=colors[0])
ax.annotate("       \n          ", (ts_rs[0], ts_zs[5]), xytext=(1.4, -0.4),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

# Langmuir probes.
lp_locs = np.array([(1.2915, -1.3655), (1.3066, -1.3655), (1.3215, -1.3655),
  (1.3365, -1.3655), (1.3515, -1.3655), (1.3665, -1.3655), (1.5003, -1.2529),
  (1.5282, -1.2529), (1.5841, -1.2529), (1.6121, -1.2529), (1.64, -1.2529),
  (1.41, -1.2529), (1.42, -1.2529)])
ax.scatter(lp_locs[:, 0], lp_locs[:, 1], s=8, color=colors[2], zorder=50)
ax.annotate("Langmuir\nProbes", lp_locs[-3], xytext=(2.2, -1.4),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

# Label the upper baffle.
ax.text(1.67, 1.10, "Upper Baffle", fontsize=fontsize, rotation=-5)

# Bold the upper divertor region and label it.
ax.plot(wall_df["rnew"][4:50], wall_df["znew"][4:50], lw=3, color="k")
ax.annotate("Upper\nDivertor", (1.25, 1.266), xytext=(1.0, 1.3),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

# Add ITF OTF labels.
ax.annotate("ITF", (2.327, -0.173), xytext=(2.419, -0.109),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize)
ax.annotate("OTF", (2.327, -0.208), xytext=(2.419, -0.30),
  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize)

# Label the outer wall.
ax.plot(wall_df["rnew"][53:59], wall_df["znew"][53:59], lw=3, color="k")
ax.text(2.22, 0.57, "Outer Wall", fontsize=fontsize, rotation=-70)

fig.tight_layout()
