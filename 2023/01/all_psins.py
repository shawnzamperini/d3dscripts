import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import numpy as np
from gadata import gadata
import MDSplus
from tqdm import tqdm
import sys

# Location of all the MP files.
root = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/"

# Option to colorcode data according to an MDSplus tag. Set to None to turn this option off.
mdsplus_tag = "DENSV2"
conn = MDSplus.Connection("atlas.gat.com")

# Isat threshold for RCP measurement. Generally something below around 0.02 might be considered noise, but there is
# no hard and fast rule here. Use at your own discretion.
isat_thresh = 0.0

data = {"psins": np.array([]), "lambdas": np.array([]), "zdata": np.array([])}
profs = {}
shots = []
for file in tqdm(os.listdir(root)):

    # Load data, filter only data above certain Isat threshold for noise.
    df = pd.read_csv(root + file, delimiter="\t")
    rcp_isat = df["Isat(A)"].values
    noise_mask = rcp_isat > isat_thresh
    rcp_r = df["R(cm)"].values[noise_mask] / 100  # cm to m
    rcp_ne = df["Ne(E18 m-3)"].values[noise_mask] * 1e18
    rcp_psin = np.square(df["Rho"].values)[noise_mask]

    # Remove any negative values.
    neg_mask = rcp_ne > 0.00
    rcp_ne = rcp_ne[neg_mask]
    rcp_r = rcp_r[neg_mask]
    rcp_lnne = np.log(rcp_ne)
    rcp_psin = rcp_psin[neg_mask]

    # Load average MDSplus tag value during plunge and store.
    shot = int(file.split("MP")[1].split("_")[0])
    shots.append(shot)
    if type(mdsplus_tag) != type(None):
        gaobj = gadata(mdsplus_tag, shot, connection=conn)
        min_time = df["Time(ms)"].min()
        max_time = df["Time(ms)"].max()
        mask = np.logical_and(gaobj.xdata >= min_time, gaobj.xdata <= max_time)
        avg_mds_tag = gaobj.zdata[mask].mean()
    else:
        avg_mds_tag = 0.0

    # Store data for radial ne profiles plot.
    profs[shot] = {"psin": rcp_psin, "ne": rcp_ne, "zdata": avg_mds_tag}

    # Do a running window of lambda_ne at each location.
    window_size = 0.03
    lambda_ne = np.zeros(len(rcp_r))
    for i in range(0, len(rcp_r)):
        window = np.full(len(rcp_r), False)
        r = rcp_r[i]
        mask = np.abs(rcp_r - r) < window_size / 2
        window[mask] = True
        z = np.polyfit(rcp_r[window], rcp_lnne[window], 1)
        p = np.poly1d(z)
        lambda_ne[i] = -1 / z[0]

    data["psins"] = np.append(data["psins"], rcp_psin)
    data["lambdas"] = np.append(data["lambdas"], lambda_ne)
    data["zdata"] = np.append(data["zdata"], np.full(len(rcp_psin), avg_mds_tag))

# Plot of the ne profiles. Need to first establish the colorbar stuff before plotting.
if type(mdsplus_tag) != type(None):
    cbar_min = sys.float_info.max; cbar_max = -sys.float_info.max
    for shot in shots:
        if profs[shot]["zdata"] < cbar_min:
            cbar_min = profs[shot]["zdata"]
        if profs[shot]["zdata"] > cbar_max:
            cbar_max = profs[shot]["zdata"]

    # Can just set the colorbar limits yourself in case something is fucking it up (i.e. DENSV2 gets goofed by 184529).
    print("Manually assigning colorbar limits!")
    cbar_min = 1e19
    cbar_max = 5e19

    # Pick your norm type.
    norm = plt.Normalize(cbar_min, cbar_max)
    #norm = colors.LogNorm(cbar_min, cbar_max)
    cmap = plt.get_cmap("inferno")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig, ax = plt.subplots(figsize=(5, 4))
for shot in shots:
    ax.plot(profs[shot]["psin"], profs[shot]["ne"], color="k", lw=2)
    if type(mdsplus_tag) != type(None):
        ax.plot(profs[shot]["psin"], profs[shot]["ne"], color=cmap(norm(profs[shot]["zdata"])), lw=1)
    else:
        ax.plot(profs[shot]["psin"], profs[shot]["ne"], color="tab:red", lw=1)
if type(mdsplus_tag) != type(None):
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(mdsplus_tag)
ax.set_xlim(1.0, 1.3)
ax.set_ylim(5e17, 2e19)
ax.set_yscale("log")
ax.grid(which="both")
ax.set_xlabel("Psin")
ax.set_ylabel("ne (m-3)")
fig.tight_layout()
fig.show()

# Plot of all the lambda_ne's.
fig, ax1 = plt.subplots(figsize=(5, 4))
if type(mdsplus_tag) != type(None):
    ax1.scatter(data["psins"], data["lambdas"] * 100, s=20, alpha=0.75, marker="^", c=data["zdata"], edgecolors="k",
                zorder=15, norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax1)
    cbar.set_label(mdsplus_tag)
else:
    ax1.scatter(data["psins"], data["lambdas"] * 100, s=20, alpha=0.75, marker="^", color="tab:red", edgecolors="k",
                zorder=15)
ax1.set_xlim([1, 1.30])
ax1.set_ylim([0, 15])
ax1.grid(zorder=5)
ax1.set_xlabel("Psin")
ax1.set_ylabel("lambda ne (cm)")
fig.tight_layout()
fig.show()
