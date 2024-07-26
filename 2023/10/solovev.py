import numpy as np
import matplotlib.pyplot as plt
import pickle
from shapely.geometry import Polygon, Point


# Following along in Mandell dissertation, Ch. 5.3, for a Solov'ev equlibrium.
R0 = 1.722
B0 = 2.002
# kappa0 = 1.695  # ellipticity (== elongation?)
kappa0 = 1.5
# qbar = 8.481  # safety factor at the axis (q1 in EFIT?)
qbar = 3.902

R = np.linspace(1.0, 2.4, 100)
Z = np.linspace(-1.5, 1.5, 100)
Rs, Zs = np.meshgrid(R, Z)
psi = B0 / (2 * R0**2 * kappa0 * qbar) * (np.square(Rs) * np.square(Zs) + kappa0**2 / 4 * np.square((np.square(Rs) - R0**2)))
psi_sep = B0 * kappa0 * R0**2 / (8 * qbar)


with open("/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle", "rb") as f:
    gfile = pickle.load(f)
efit_Rs, efit_Zs = np.meshgrid(gfile["R"], gfile["Z"])
efit_psi = gfile["PSIRZ"]
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
mask = np.full(Rs.shape, False)
efit_mask = np.full(efit_Rs.shape, False)
for i in range(0, Rs.shape[0]):
    for j in range(0, Rs.shape[1]):
        mask[i, j] = wall.contains(Point(Rs[i, j], Zs[i, j]))
for i in range(0, efit_Rs.shape[0]):
    for j in range(0, efit_Rs.shape[1]):
        efit_mask[i, j] = wall.contains(Point(efit_Rs[i, j], efit_Zs[i, j]))
psi_masked = np.ma.masked_where(~mask, psi)
efit_psi_masked = np.ma.masked_where(~efit_mask, efit_psi)

levels = np.linspace(gfile["SIBRY"], 0.35, 5)
fig, ax = plt.subplots()
ax.contour(efit_Rs, efit_Zs, efit_psi_masked, colors="k", levels=levels)
ax.contour(efit_Rs, efit_Zs, efit_psi_masked, colors="k", levels=[gfile["SIBRY"]], linewidths=3)
ax.contour(Rs, Zs, psi_masked, colors="r", levels=levels)
ax.contour(Rs, Zs, psi_masked, colors="r", levels=[psi_sep], linewidths=3)
ax.plot(rwall, zwall, color="k", lw=2)
ax.set_ylim(-1.55, 1.55)
ax.set_xlim(0.85, 2.5)
ax.plot([-99, -99], [-99, -99], color="k", label="EFIT")
ax.plot([-99, -99], [-99, -99], color="r", label="Solov'ev")
ax.legend(loc="upper right")
ax.set_aspect("equal")
fig.tight_layout()
fig.show()
