import pickle
import matplotlib.pyplot as plt
import matplotlib.path as mplpath
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.patches as patches


shot = 167247

if shot == 167247:
    nearsol_end = 1.2
    farsol_end  = 8.5
elif shot == 167277:
    nearsol_end = -0.5
    farsol_end  = 2.3

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# TS data for later.
ts_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/setup-files/ts167247_final_v2.xlsx"
ts_df = pd.read_excel(ts_path)

root = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/{}/setup-files/{}_".format(shot, shot)

def load_pickle(path):
    with open(path, "rb") as f:
        var = pickle.load(f)
    return var

psirz = load_pickle(root+"psirz")
rbbbs = load_pickle(root+"rbbbs")
zbbbs = load_pickle(root+"zbbbs")
rlim  = load_pickle(root+"rlim")
zlim  = load_pickle(root+"zlim")
rdim  = load_pickle(root+"rdim")
zdim  = load_pickle(root+"zdim")
rleft = load_pickle(root+"rleft")
zmid  = load_pickle(root+"zmid")
sibry = load_pickle(root+"sibry")

# Fix outer wall jankyness.
rlim_old = rlim.copy()
zlim_old = zlim.copy()
m = (zlim[56] - zlim[55]) / (rlim[56] - rlim[55])
r_fix = rlim[66]
z_fix = m * (rlim[66] - rlim[55]) + zlim[55]
#rlim[56:66] = np.full(len(rlim[56:66]), rlim[66])
#zlim[56:66] = np.full(len(zlim[56:66]), z_fix)
rlim = np.delete(rlim, np.arange(57, 66, dtype=np.int))
zlim = np.delete(zlim, np.arange(57, 66, dtype=np.int))
rlim = np.insert(rlim, 57, r_fix)
zlim = np.insert(zlim, 57, z_fix)

R, Z = np.meshgrid(np.linspace(0, rdim, psirz.shape[0])+rleft, np.linspace(0, zdim, psirz.shape[0])-zdim/2)

# Mask the points outside of the vessel boundary. This command is bangin'.
bbpath = mplpath.Path(list(zip(rlim, zlim)))
bbpath_mask = ~bbpath.contains_points(np.array(list(zip(R.flatten(), Z.flatten()))))
psirz_masked = np.ma.masked_array(psirz.flatten(), mask=bbpath_mask).reshape(psirz.shape)

psin = -psirz_masked / sibry
#nearsol_end = 1.2
#farsol_end  = 8
#core    = np.ma.masked_array(psin, mask=np.logical_and(psin>-1, Z>-1.103))
core    = np.ma.masked_array(psin, mask=psin>-1)
sol     = np.ma.masked_array(psin, mask=psin<-1)
nearsol = np.ma.masked_array(psin, mask=~np.logical_and(psin>-1, psin<nearsol_end+1))
farsol  = np.ma.masked_array(psin, mask=~np.logical_and(psin>=nearsol_end, psin<farsol_end+0.2))
wallsol = np.ma.masked_array(psin, mask=~np.logical_and(psin>=farsol_end, R>=1.6))

# Extra lines for the legs.
xp = zbbbs.argmin()
#leg1_m = (zbbbs[xp] - zbbbs[xp-1]) / (rbbbs[xp] - rbbbs[xp-1])
#leg2_m = (zbbbs[xp] - zbbbs[xp+1]) / (rbbbs[xp] - rbbbs[xp+1])
leg1_r = np.linspace(rbbbs[xp], 1.300,  100)
leg1_z = np.linspace(zbbbs[xp], -1.36,  100)
leg2_r = np.linspace(rbbbs[xp], 1.015,  100)
leg2_z = np.linspace(zbbbs[xp], -1.113, 100)

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14

fig, ax = plt.subplots(figsize=(4.34, 7.19))
cont0 = ax.contour(R, Z, core, linestyles="solid", colors="green")
cont0 = ax.contour(R, Z, sol, colors="blue")
#cont1 = ax.contourf(R, Z, nearsol, cmap=mpl.colors.ListedColormap([colors[2]]))
#cont3 = ax.contourf(R, Z, farsol,  cmap=mpl.colors.ListedColormap([colors[3]]))
#cont4 = ax.contourf(R, Z, wallsol, cmap=mpl.colors.ListedColormap([colors[4]]))
#cbar = fig.colorbar(cont)
ax.plot(rbbbs, zbbbs, color="k", lw=4)
ax.plot(leg1_r, leg1_z, color="k", lw=4)
ax.plot(leg2_r, leg2_z, color="k", lw=4)
ax.plot(rlim, zlim, color="k", lw=2)
ax.set_aspect("equal")
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="both", which="both", labelleft=False, labelbottom=False, left=False, bottom=False)

# Add CP.
#cp_r_tip = 2.25
#cp_z = -0.188
#cp_width = 0.03
#rect_anchor = (cp_r_tip, cp_z - cp_width / 2)  # Z location of anchor just the bottom of CP
#rect = patches.Rectangle(rect_anchor, 2.377-cp_r_tip, cp_width, ec="k", color=colors[2])
#rect = patches.Rectangle(rect_anchor, 2.3755-cp_r_tip, cp_width, ec="k", color=colors[2])
#ax.add_patch(rect)

# Coordinates for Thomson scattering. Get them from the TS file.
#ts_core = ts_df[ts_df["System"] == "core"]
#ts_zs = ts_core["Z (m)"].unique()
#ts_rs = np.full(len(ts_zs), ts_core["R (m)"].unique()[0])
#ax.scatter(ts_rs, ts_zs, s=3, color=colors[0])
#ax.annotate("Thomson\nScattering", (ts_rs[0], ts_zs[35]), xytext=(1.8, -0.4),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

#ts_div = ts_df[ts_df["System"] == "divertor"]
#ts_zs = ts_div["Z (m)"].unique()
#ts_rs = np.full(len(ts_zs), ts_div["R (m)"].unique()[0])
#ax.scatter(ts_rs, ts_zs, s=3, color=colors[0])
#ax.annotate("       \n          ", (ts_rs[0], ts_zs[5]), xytext=(1.8, -0.4),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

# Langmuir probes.
#lp_locs = np.array([(1.2915, -1.3655), (1.3066, -1.3655), (1.3215, -1.3655),
#  (1.3365, -1.3655), (1.3515, -1.3655), (1.3665, -1.3655), (1.5003, -1.2529),
#  (1.5282, -1.2529), (1.5841, -1.2529), (1.6121, -1.2529), (1.64, -1.2529),
#  (1.41, -1.2529), (1.42, -1.2529)])
#ax.scatter(lp_locs[:, 0], lp_locs[:, 1], s=8, color="k", zorder=50)
#ax.annotate("Langmuir\nProbes", lp_locs[-3], xytext=(2.2, -1.4),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

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
                                   edgecolor=edgecolor, zorder=99)
shelf_rect  = patches.Rectangle(shelf_xy, width = shelf_width,
                                   height=tile_height,
                                   facecolor=facecolor,
                                   edgecolor=edgecolor, zorder=99)
#ax.add_patch(floor_rect)
#ax.add_patch(shelf_rect)

# Add metal ring annotation.
#ax.annotate("Floor Ring", (floor_xy[0] + floor_width / 2, floor_xy[1]), xytext=(1.2, -1.55),
#  arrowprops=dict(facecolor="black", arrowstyle="-"),
#  fontsize=fontsize)
#ax.annotate("Shelf Ring", (shelf_xy[0] + shelf_width / 2, shelf_xy[1]), xytext=(1.55, -1.45),
#  arrowprops=dict(facecolor="black", arrowstyle="-"),
#  fontsize=fontsize)

# Label the upper baffle.
#ax.text(1.64, 1.10, "Upper Baffle", fontsize=fontsize, rotation=-5)

# Bold the upper divertor region and label it.
#ax.plot(rlim[4:49], zlim[4:49], lw=3, color="k")
#ax.annotate("Upper\nDivertor", (1.25, 1.266), xytext=(1.0, 1.3),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize, ha="center")

# Add ITF OTF labels.
#ax.annotate("ITF", (2.327, -0.173), xytext=(2.419, -0.109),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize)
#ax.annotate("OTF", (2.327, -0.208), xytext=(2.419, -0.30),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize)

# Label the outer wall.
#ax.plot(rlim[54:62], zlim[54:62], lw=3, color="k")
#ax.text(2.19, 0.57, "Outer Wall", fontsize=fontsize, rotation=-70)
#ax.text(0.9, -0.4, "Inner Wall", fontsize=fontsize, rotation=90)

# Annotate the SOL regions.
#ax.annotate("Near-\nSOL", (1.323, -1.14), xytext=(1.244, -0.75),
#  arrowprops=dict(facecolor="black", arrowstyle="-"), fontsize=fontsize)

#ax.text(1.2, -0.634, "Near-SOL", fontsize=fontsize, rotation=-70)
#ax.text(1.12, 0.87, "Far-SOL", fontsize=fontsize, rotation=35)
#ax.text(1.9, 0.7, "Wall-SOL", fontsize=fontsize, rotation=-50)

fig.tight_layout()
fig.show()
