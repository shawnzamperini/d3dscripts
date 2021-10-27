import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm
from matplotlib import ticker, cm
from matplotlib.colors import LogNorm


nR = 100
nZ = 100
col_names = ["R[m]", "Z[m]", "N_toroidal", "connection length [km]",
  "psimin (penetration depth)", "psimax", "psiav", "FL pitch angle",
  "FL yaw angle", "theta", "psi"]
root = "/Users/zamperini/Documents/d3d_work/mafot_187122/"

def read_df(fname):
    df = pd.read_csv(root+fname, delimiter="\t", skiprows=52, names=col_names)
    return df

def reshape(df, col):
    return np.reshape(df[col].values, (nR, nZ))

# Some variables are the same for every run, like R and Z, so extract those
# ahead of time since we don't have to copy those over each time.
df0 = read_df("lam_torang0.dat")
full_df = pd.DataFrame()
full_df["R"] = df0["R[m]"]
full_df["Z"] = df0["Z[m]"]

for tor in tqdm(range(0, 360)):

    df = read_df("lam_torang{}.dat".format(tor))

    full_df["{}_conn".format(tor)] = df["connection length [km]"] * 1000
    full_df["{}_psin".format(tor)] = df["psi"]


R = reshape(full_df, "R")
Z = reshape(full_df, "Z")
conn = reshape(full_df, "0_conn")
psin = reshape(full_df, "0_psin")
mask = psin < 1
R = np.ma.masked_array(R, mask=mask)
Z = np.ma.masked_array(Z, mask=mask)
conn = np.ma.masked_array(conn, mask=mask)
psin = np.ma.masked_array(psin, mask=mask)

fig, ax = plt.subplots(figsize=(3,6))
ax.set_facecolor("grey")
c = ax.contourf(R, Z, conn, norm=LogNorm(), levels=np.geomspace(0.1, 100, 10), cmap="magma")
cbar = fig.colorbar(c)
ax.set_aspect("equal")
fig.show()
