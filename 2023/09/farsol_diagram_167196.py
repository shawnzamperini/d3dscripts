import pickle
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, Point


with open("/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle", "rb") as f:
    gfile = pickle.load(f)

Rs, Zs = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]
rwall = list(gfile["RLIM"])
zwall = list(gfile["ZLIM"])

# Modify wall to remove the helicon thing using point-slope formula.
x1 = rwall[55]
x2 = rwall[56]
y1 = zwall[55]
y2 = zwall[56]
m = (y2 - y1) / (x2 - x1)
newr = rwall[66]
newz = m * (newr - x1) + y1

# Remove bad points, plug in good ones.
del(rwall[56:66])
del(zwall[56:66])
rwall.insert(56, newr)
zwall.insert(56, newz)

wall = Polygon(zip(rwall, zwall))
mask1 = np.full(Rs.shape, False)
mask2 = np.full(Rs.shape, False)
mask3 = np.full(Rs.shape, False)
for i in range(0, Rs.shape[0]):
    for j in range(0, Rs.shape[1]):
        if Rs[i, j] > 1.625:
            mask1[i, j] = wall.contains(Point(Rs[i, j], Zs[i, j]))
        mask2[i, j] = wall.contains(Point(Rs[i, j], Zs[i, j]))
        mask3[i, j] = np.logical_and(wall.contains(Point(Rs[i, j], Zs[i, j])), Zs[i, j] > gfile["Zx1"])
psin_masked1 = np.ma.masked_where(~mask1, psin)
psin_masked2 = np.ma.masked_where(~mask2, psin)
psin_masked3 = np.ma.masked_where(~mask3, psin)

lp_locs = np.array([[1.5003, -1.2529],
                    [1.5282, -1.2529],
                    [1.5841, -1.2529],
                    [1.6121, -1.2529],
                    [1.64, -1.2529],
                    [1.6679, -1.2529],
                    [1.41, -1.2529],
                    [1.42, -1.2529],
                    [1.078, -1.278]])

psin_window = 1.127
fontsize = 14
alpha = 0.5
fig, ax = plt.subplots(figsize=(5, 6))
ax.plot(rwall, zwall, color="k")
ax.contourf(Rs, Zs, psin_masked2, levels=[0, 99], colors="tab:pink", alpha=alpha)
ax.contourf(Rs, Zs, psin_masked3, levels=[0.0, 1.0], colors="w")
ax.contourf(Rs, Zs, psin_masked3, levels=[0.0, 1.0], colors="tab:red", alpha=alpha)
ax.contourf(Rs, Zs, psin_masked1, levels=[psin_window, 1.8], colors="white")
ax.contourf(Rs, Zs, psin_masked1, levels=[psin_window, 1.8], colors="tab:green", alpha=alpha)
ax.contour(Rs, Zs, psin_masked1, levels=[psin_window], colors="k", linewidths=2)
ax.contour(Rs, Zs, psin_masked2, levels=[1.00], colors="k", linewidths=3)
ax.plot([1.94, 1.94], [-0.068, 0.827], color="k")
ax.annotate("Thomson\nScattering", (1.94, 0.4), xytext=(1.57, 0.4), color="k",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="k", lw=1), fontsize=fontsize, ha="center")
ax.plot([2.230, 2.366], [-0.185, -0.185], color="k")
ax.annotate("RCP", (2.28, -0.185), xytext=(1.85, -0.35), color="k",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="k"), fontsize=fontsize, ha="center")
ax.scatter(lp_locs[:, 0], lp_locs[:, 1], s=8, color="k", zorder=50)
ax.annotate("Langmuir\nProbes", lp_locs[2], xytext=(2.2, -1.45), color="k",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="k"), fontsize=fontsize, ha="center")
ax.plot([1.404, 1.454], [-1.25, -1.25], color="tab:cyan", lw=3, zorder=55)
ax.plot([1.32, 1.37], [-1.363, -1.363], color="tab:cyan", lw=3, zorder=55)
ax.annotate("W Rings", (1.425, -1.25), xytext=(1.6, -1.55), color="tab:cyan",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="tab:cyan"), fontsize=fontsize, ha="center")
ax.annotate("W Rings", (1.345, -1.363), xytext=(1.6, -1.55), color="tab:cyan",
            arrowprops=dict(facecolor="black", arrowstyle="-", edgecolor="tab:cyan"), fontsize=fontsize, ha="center")
ax.set_aspect("equal")
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.text(1.23, -0.35, "TGYRO", fontsize=fontsize, color="tab:red", fontweight="semibold")
ax.text(1.06, 0.95, "DIVIMP +\nOSM-EIRENE", fontsize=12, color="tab:pink", fontweight="semibold")
ax.text(1.85, 0.70, "3DLIM", fontsize=fontsize, rotation=-45, color="tab:green", fontweight="semibold")
ax.text(1.67, 1.08, "#167196", fontsize=10, rotation=-5)
ax.set_xticks([])
ax.set_yticks([])
fig.tight_layout()
fig.show()
