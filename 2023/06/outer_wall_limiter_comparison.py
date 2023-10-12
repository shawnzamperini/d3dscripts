from LimWallToolkit import LimWallToolkit
import matplotlib.pyplot as plt
import pickle
import numpy as np


lwt = LimWallToolkit()
wall1 = lwt.read_3d_wall("/Users/zamperini/Documents/d3d_work/mafot_files/wall_files/mafot_wall_gap_torlim_v2.dat")
wall2 = lwt.read_3d_wall("/Users/zamperini/Documents/d3d_work/mafot_files/wall_files/mafot_3d_wall.dat")

# Get gfile for plotting flux surfaces.
with open("/Users/zamperini/Documents/d3d_work/mafot_files/186257/186257_3500.pickle", "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]


fig, ax = plt.subplots()
ax.contour(Rs, Zs, psin, colors="k", alpha=0.10)
ax.contour(Rs, Zs, psin, colors="k", levels=[1.15])
ax.plot(wall1[130][0], wall1[130][1], color="g", label="New", lw=3)
ax.plot(wall2[130][0], wall2[130][1], color="r", linestyle="-", label="Current", lw=3)
ax.set_aspect("equal")
ax.set_xlabel("R (m)", fontsize=12)
ax.set_ylabel("Z (m)", fontsize=12)
ax.set_title("#186257 (SVR)", fontsize=12)
ax.set_ylim([-0.5, 0.5])
ax.set_xlim([2.1, 2.5])
#ax.legend()
fig.tight_layout()
fig.show()