import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata, interp2d
from tqdm import tqdm
import pandas as pd


shot = 171646  # Either 171650 or 171646
path = "/Users/zamperini/Documents/d3d_work/ts_{}.pickle".format(shot)
with open(path, "rb") as f:
    ts = pickle.load(f)

time = 2500
trange = 250
rsep = 1.94
zsep = 0.473
# idx = np.argmin(np.abs(ts["core"]["time"] - time))
idx = np.logical_and(ts["core"]["time"] > (time - trange), ts["core"]["time"] < time + trange)
psin = ts["core"]["psin"][:, idx].flatten()
te = ts["core"]["te"][:, idx].flatten()
ne = ts["core"]["ne"][:, idx].flatten()
r = np.tile(ts["core"]["r"], idx.sum())
z = np.tile(ts["core"]["z"], idx.sum())

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
gfile_path = "/Users/zamperini/Documents/d3d_work/{}_2500.pickle".format(shot)
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psins = gfile["PSIRZ_NORM"]
Rtrunc = Rs > gfile["RMAXIS"]
coords = tuple(zip(psin, np.full(len(psin), gfile["ZMAXIS"])))
Romps = []
for i in tqdm(range(0, len(coords))):
    x = coords[i][0]
    y = coords[i][1]
    Romps.append(float(griddata((psins[Rtrunc].flatten(), Zs[Rtrunc].flatten()), Rs[Rtrunc].flatten(), (x, y))))
# R_omps = griddata((psins[Rtrunc].flatten(), Zs[Rtrunc].flatten()), Rs[Rtrunc].flatten(), coords)
# f_Romp = interp2d(psins[Rtrunc], Zs[Rtrunc], Rs[Rtrunc])
Rsep_omp = float(griddata((psins[Rtrunc].flatten(), Zs[Rtrunc].flatten()), Rs[Rtrunc].flatten(), (1.0, gfile["ZMAXIS"])))
rmrsomp = np.array(Romps) - Rsep_omp


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
ax1.axvline(0.0, color="k", linestyle="--")
ax2.axvline(0.0, color="k", linestyle="--")
# ax1.plot(psin, te, marker=".")
# ax2.plot(psin, ne, marker=".")
ax1.scatter(rmrsomp, te, alpha=0.4)
ax2.scatter(rmrsomp, ne, alpha=0.4)
ax1.set_xlim([-0.05, 0.15])
ax1.set_ylim([0, 200])
ax1.set_xlabel("R-Rsep @ OMP (m)")
ax2.set_xlabel("R-Rsep @ OMP (m)")
ax1.set_ylabel("Te (eV)")
ax2.set_ylabel("ne (m-3)")
# fig.suptitle("171650 {} +/- {}".format(time, trange))
fig.tight_layout()
fig.show()

# Save as csv.
df = pd.DataFrame({"psin":psin, "R @ OMP (m)":Romps, "R-Rsep @ OMP (m)":rmrsomp, "Te (eV)":te, "ne (m-3)":ne})
df.to_csv("ts{}_{}-{}.csv".format(shot, time-trange, time+trange))