import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplpath
from LimWallToolkit import LimWallToolkit


gfile_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/gfile_data_for_rcp/184527_3400"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)

lwt = LimWallToolkit()
wall_path = "/Users/zamperini/Documents/d3d_work/184527/mafot/mafot_3D_wall.dat"
wall = lwt.read_3d_wall(wall_path)
wall_coords = wall[243]

# Fix issue in inner wall coordinates.
del wall_coords[0][-3:]
del wall_coords[1][-3:]
del wall_coords[0][:10]
del wall_coords[1][:10]
wall_coords[0].append(wall_coords[0][0])
wall_coords[1].append(wall_coords[1][0])

fig, ax = plt.subplots(figsize=(4,6))

ax.plot(wall_coords[0], wall_coords[1], color="k", zorder=3)
#ax.scatter(wall_coords[0][8], wall_coords[1][8], color="r")

# Easiest way here would be to plot a normal contour plot at
# specified psin steps, and to mask anything outside of the vessel.
R, Z = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["Psin"]
bbpath = mplpath.Path(list(zip(wall_coords[0], wall_coords[1])))
bbpath_mask = ~bbpath.contains_points(np.array(list(zip(R.flatten(),
  Z.flatten()))))
psin_masked = np.ma.masked_array(psin.flatten(),
  mask=bbpath_mask).reshape(psin.shape)
ax.contour(gfile["R"], gfile["Z"], psin_masked, colors="k",
  levels=np.geomspace(1.0, 1.15, 7), zorder=2, linewidths=1)
ax.contour(gfile["R"], gfile["Z"], psin_masked, colors="k",
  levels=[1.0], linewidths=2, zorder=1)

ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", which="both", bottom=False, top=False,
  left=False, right=False, labelbottom=False, labelleft=False)
ax.set_aspect("equal")
fig.show()
