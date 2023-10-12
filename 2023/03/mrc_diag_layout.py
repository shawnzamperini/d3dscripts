# Plot the layout of the MRC diagnostics that are used for OEDGE constraining for 167196.
import matplotlib.pyplot as plt
import pickle
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# Load flux data from gfile.
gfile_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/167195_4460.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
wallr = gfile["RLIM"]
wallz = gfile["ZLIM"]

# Create a mask so we don't plot flux surfaces outside the vessel wall.
wall_poly = Polygon(zip(wallr, wallz))
mask = np.full(Rs.shape, True)
for i in range(0, mask.shape[0]):
    for j in range(0, mask.shape[1]):
        p = Point(Rs[i, j], Zs[i, j])
        if wall_poly.contains(p):
            mask[i, j] = False
Rs = np.ma.masked_where(mask, Rs)
Zs = np.ma.masked_where(mask, Zs)
psin = np.ma.masked_where(mask, psin)

fig, ax = plt.subplots(figsize=(3, 5))
ax.set_aspect("equal")
ax.plot(wallr, wallz, color="k")
ax.contour(Rs, Zs, psin, levels=[1.0], colors="k")
ax.contour(Rs, Zs, psin, levels=[1.05, 1.1, 1.15, 1.2, 1.25, 1.3], alpha=0.2, colors="k")

ax.plot([1.94, 1.94], [-0.068, 0.827], color="k")
ax.annotate("Thomson\nScattering", (1.94, 0.4), xytext=(1.57, 0.4), color="k",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="k", lw=0), fontsize=12, ha="center")

ax.plot([2.230, 2.366], [-0.185, -0.185], color="k")
ax.annotate("RCP", (2.28, -0.185), xytext=(1.85, -0.35), color="k",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="k"), fontsize=12, ha="center")

lp_locs = np.array([[1.5003, -1.2529],
                    [1.5282, -1.2529],
                    [1.5841, -1.2529],
                    [1.6121, -1.2529],
                    [1.64, -1.2529],
                    [1.6679, -1.2529],
                    [1.41, -1.2529],
                    [1.42, -1.2529],
                    [1.078, -1.278]])
ax.scatter(lp_locs[:, 0], lp_locs[:, 1], s=8, color="k", zorder=50)
ax.annotate("Langmuir\nProbes", lp_locs[2], xytext=(2.2, -1.45), color="k",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="k"), fontsize=12, ha="center")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.tick_params(axis="both", which="both", labelleft=False, labelbottom=False, left=False, bottom=False)
fig.tight_layout()
fig.show()
