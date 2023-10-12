import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append("../../2022/12")
import BlobbyFarSOL

shot = 190484
plunge = 2
window_size = 0.015

bfs = BlobbyFarSOL.BlobbyFarSOL(shot, plunge)
bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190484_2.tab", 20, 2)
r = bfs.rcp_r / 100  # cm to m
ne = bfs.rcp_ne
lnne = np.log(ne)

ex_plots = {}
for i in range(0, len(r)):
    window = np.full(len(r), False)
    mask = np.abs(r - r[i]) < window_size / 2
    window[mask] = True
    z = np.polyfit(r[window], lnne[window], 1)
    p = np.poly1d(z)
    lambda_ne = -1 / z[0]

    # Save a few examples.
    if i in [4, 9, 14]:
        ex_plots[i] = {"ri":r[i], "nei":np.exp(lnne[i]), "r":r[window], "ne":np.exp(lnne[window]), "lambda_ne":lambda_ne}

colors = ["tab:red", "tab:cyan", "tab:pink"]
count = 0
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.scatter(r, ne, color="k")
for k, v in ex_plots.items():
    ax1.scatter(v["r"], v["ne"], color=colors[count], edgecolors="k")
    ax1.scatter(v["ri"], v["nei"], color=colors[count], edgecolors="k", marker="*", s=200)
    ax1.plot(v["r"], v["nei"] * np.exp(-(v["r"]-v["ri"]) / v["lambda_ne"]), color=colors[count])
    ax1.annotate(r"$\mathdefault{\lambda_{n}}$" + " = {:.2f} cm".format(v["lambda_ne"] * 100), (v["ri"], v["nei"]),
                 color=colors[count], textcoords="offset points", xytext=(10, 13), fontsize=10, arrowprops={"width":2, "headwidth":0, "color":colors[count]})
    count += 1
ax1.set_yscale("log")
ax1.grid(alpha=0.25, which="both")
#ax1.set_yticks([2e18, 4e18, 6e18, 8e18])
#ax1.set_yticklabels([])
#ax1.set_yticklabels(["{}".format(i) + r"$\mathdefault{ \times 10^{18}}$" for i in [2, 4, 6, 8]])
ax1.set_ylabel(r"$\mathdefault{n_e\ (m^{-3})}$", fontsize=12)
ax1.set_xlabel("R (m)", fontsize=12)
ax1.text(2.277, 7e18, "#190484", fontsize=12, bbox={"facecolor":"w"})
fig.tight_layout()
fig.show()