# This script will generate a long plot just showing the exponential and flat
# impurity distributions for the W source.
import numpy as np
import matplotlib.pyplot as plt


xmin = 0
xmax = 6.99 + 1.6
lamb = 0.65
x = np.linspace(xmin, xmax, 100)
A = lamb * (np.exp(xmax/lamb) - np.exp(xmin/lamb))
y = np.exp(x/lamb) / A

A2 = (xmax - xmin)
y2 = np.full(len(x), 1/A2)
y2 = np.full(len(x), 0.5)

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fontsize = 14
lw = 5

fig, ax = plt.subplots(figsize=(6, 2))
ax.plot(x, y, label="Unfavorable", lw=lw, color=colors[2])
#ax.plot(x, y2, label="Favorable", lw=lw, color=colors[3])
ax.set_xlabel("Distance from outer target", fontsize=fontsize)
#ax.tick_params(axis="both", which="both", bottom=False, left=False, labelbottom=False, labelleft=False)
#ax.legend(fontsize=fontsize, loc="upper center")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
fig.tight_layout()
fig.show()
