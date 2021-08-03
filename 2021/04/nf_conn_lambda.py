import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

lambdas = [9.6664,6.7665,4.7501,3.5282,3.3296,3.5527,4.5319,4.6239,4.7675]
conns = [68.42,62.38,57.53,53.64,48.05,45.68,44.58,44.5,45.82]

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 6))
fontsize = 14
lw = 5

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(conns, lambdas, color=colors[4], s=80, edgecolors="k", linewidths=1)
ax.grid()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlabel("Connection Length (m)", fontsize=fontsize)
ax.set_ylabel("Exponential Falloff (m)", fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=12)
fig.tight_layout()
fig.show()
