import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata, interp2d
from tqdm import tqdm
import pandas as pd
from freeqdsk import geqdsk
import pgkylUtil_adios2 as pgu


shot = 163150
path = "/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-163150-v1/ts_{}.pickle".format(shot)
with open(path, "rb") as f:
    ts = pickle.load(f)

time = 1500
trange = 100
# rsep = 1.94
# zsep = 0.473
# idx = np.argmin(np.abs(ts["core"]["time"] - time))
idx = np.logical_and(ts["core"]["time"] > (time - trange), ts["core"]["time"] < time + trange)
psin = ts["core"]["psin"][:, idx].flatten()
te = ts["core"]["te"][:, idx].flatten()
ne = ts["core"]["ne"][:, idx].flatten()
# r = np.tile(ts["core"]["r"], idx.sum())
# z = np.tile(ts["core"]["z"], idx.sum())
r = np.repeat(ts["core"]["r"], idx.sum())
z = np.repeat(ts["core"]["z"], idx.sum())


# Not interested in anything where psin = 0.
mask = psin != 0
psin = psin[mask]
te = te[mask]
ne = ne[mask]
r = r[mask]
z = z[mask]

sort_idx = np.argsort(psin)
psin = psin[sort_idx]
te = te[sort_idx]
ne = ne[sort_idx]
r = r[sort_idx]
z = z[sort_idx]

# Will want to map to the outboard midplane.
gfile_path = "/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-163150-v1/g163150.01500"
with open(gfile_path, "r") as f:
    gfile = geqdsk.read(f)
gfile_R = np.linspace(gfile["rleft"], gfile["rleft"]+gfile["rdim"], gfile["nx"])
gfile_Z = np.linspace(gfile["zmid"] - gfile["zdim"]/2, gfile["zmid"] + gfile["zdim"]/2, gfile["ny"])
gfile_Rs, gfile_Zs = np.meshgrid(gfile_R, gfile_Z, indexing="ij")
gfile_psis = gfile["psi"]
rlim = gfile["rlim"]
zlim = gfile["zlim"]

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.contourf(gfile_Rs, gfile_Zs, gfile_psis)
ax.contour(gfile_Rs, gfile_Zs, gfile_psis, levels=[gfile["sibdry"]], colors="r")
ax.plot(rlim, zlim, color="k")
fig.tight_layout()
fig.show()

# -----------
# Need the Gkeyll geometry. Copied form Tess's script plot_equil_posD.py.
Raxis   = 0.86
gridNodes = pgu.getRawData('/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-163150-v1/d3d-shot163150-axisym_grid.bp')
Rp, Zp, zp = list(), list(), list()

simDim = np.size(np.shape(gridNodes)) - 1
slices = [np.s_[0:1]] * (simDim + 1)
slices[simDim] = 0
for d in range(simDim):
    slices[d] = np.s_[0:]
numCells = np.shape(gridNodes[tuple(slices)])

nodeList = np.reshape(gridNodes, [np.product(numCells), 3])

# .Compute field line following coordinate z to color the nodes.
# .Below we refer to Cartesian coordinates (X,Y,Z), cylindrical
# .coordinates (R,phi,Z) and toroidal coordinates (r,theta,phi).
# .The field aligned coordinate is z.
X, Y, gkyl_Z = nodeList[:, 0], nodeList[:, 1], nodeList[:, 2]
phi = np.arctan2(Y, X)
gkyl_R = X / np.cos(phi)
theta = np.arctan2(gkyl_Z, (gkyl_R - Raxis))
gkyl_r = (gkyl_R - Raxis) / np.cos(theta)  # = R*cos(phi)
eps = gkyl_r / Raxis

Rp.append(gkyl_R)
Zp.append(gkyl_Z)

# Shawn additions. Find the psin of every gkyl grid point.
gkyl_psi = griddata((gfile_Rs.flatten(), gfile_Zs.flatten()), gfile_psis.flatten(), (gkyl_R.flatten(), gkyl_Z.flatten()))

# -----------

ts_psi = griddata((gfile_Rs.flatten(), gfile_Zs.flatten()), gfile_psis.flatten(), (r, z))
Rtrunc = gfile_Rs > gfile["rmagx"]
coords = tuple(zip(ts_psi, np.full(len(ts_psi), gfile["zmagx"])))
Romps = []
for i in tqdm(range(0, len(coords))):
    x = coords[i][0]
    y = coords[i][1]
    Romps.append(float(griddata((gfile_psis[Rtrunc].flatten(), gfile_Zs[Rtrunc].flatten()), gfile_Rs[Rtrunc].flatten(), (x, y))))
# R_omps = griddata((psins[Rtrunc].flatten(), Zs[Rtrunc].flatten()), Rs[Rtrunc].flatten(), coords)
# f_Romp = interp2d(psins[Rtrunc], Zs[Rtrunc], Rs[Rtrunc])
Rsep_omp = float(griddata((gfile_psis[Rtrunc].flatten(), gfile_Zs[Rtrunc].flatten()), gfile_Rs[Rtrunc].flatten(), (gfile["sibdry"], gfile["zmagx"])))
rmrsomp = np.array(Romps) - Rsep_omp


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
ax1.axvspan(gkyl_psi.min(), gkyl_psi.max(), color="tab:purple", alpha=0.3)
ax2.axvspan(gkyl_psi.min(), gkyl_psi.max(), color="tab:purple", alpha=0.3)
ax1.axvline(0.0, color="k", linestyle="--")
ax2.axvline(0.0, color="k", linestyle="--")
ax1.scatter(ts_psi, te, alpha=0.4, color="tab:red", edgecolors="k")
ax2.scatter(ts_psi, ne, alpha=0.4, color="tab:red", edgecolors="k")
ax1.set_xlim([-0.05, None])
ax1.set_ylim([0, 800])
ax2.set_ylim([0.5e19, 3e19])
ax1.text(0.45, 0.7, "Gkeyll domain", transform=ax1.transAxes, color="tab:purple")
ax1.set_xlabel(r"$\mathdefault{\psi}$ (W/rad)")
ax2.set_xlabel(r"$\mathdefault{\psi}$ (W/rad)")
ax1.set_ylabel("Te (eV)")
ax2.set_ylabel("ne (m-3)")
fig.tight_layout()
fig.show()

# Save as csv.
df = pd.DataFrame({"psin":psin, "R @ OMP (m)":Romps, "R-Rsep @ OMP (m)":rmrsomp, "Te (eV)":te, "ne (m-3)":ne})
df.to_csv("ts{}_{}-{}.csv".format(shot, time-trange, time+trange))