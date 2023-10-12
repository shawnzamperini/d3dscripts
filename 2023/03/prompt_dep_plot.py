import matplotlib.pyplot as plt
import numpy as np


a_gut = 1.485
b_gut = -0.56
a_zmp = 0.9
b_zmp = -0.4

x = np.geomspace(0.01, 5, 100)
gut = np.exp(-a_gut * np.power(x, b_gut))
zmp = np.exp(-a_zmp * np.power(x, b_zmp))

ring_color = "firebrick"
zmp_color = "seagreen"
gut_color = "royalblue"

fig, ax = plt.subplots(figsize=(5, 4))

ax.axvspan(2e-2, 5e-2, color=ring_color, alpha=0.4, zorder=5, edgecolor=None)
ax.plot(x, zmp, label="Zamperini", color=zmp_color, lw=3, zorder=20)
ax.plot(x, gut, label="Guterl", color=gut_color, lw=3, zorder=15)

ax.legend(fontsize=12)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel(r"$\mathdefault{1-f_{redep}}$", fontsize=14)
ax.set_xlabel(r"$\mathdefault{x = \lambda_{iz} / \lambda_{sheath}}$", fontsize=14)
ax.set_xlim(0.01, 5)
ax.set_ylim([1e-4, 1])
ax.grid(alpha=0.3, which="both", zorder=10)
ax.text(0.027, 0.1, "W Ring", fontsize=12, zorder=30, color=ring_color, rotation=-90)  # bbox={"facecolor":"tab:red"},

fig.tight_layout()
fig.show()

