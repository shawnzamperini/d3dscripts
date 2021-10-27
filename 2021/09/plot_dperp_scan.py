# Quick script to plot a scan in the diffusion coefficient for a couple of
# methane shots. I used the same naming convention for each (d3d-{shot}-inj-00{2-6})
# so either can work.
import sys
sys.path.append("/Users/zamperini/github/utk-fusion/oedge/")
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


shot = 184267  # Either 184267 or 184527.
root = "/Users/zamperini/Documents/d3d_work/{}/".format(shot)

c_dens = {}
for i in range(2, 7):
    path = root + "d3d-{}-inj-00{}.nc".format(shot, i)
    op = oedge_plots.OedgePlots(path)
    x, y = op.along_ring(60, "DDLIMS", charge="all", plot_it=False)
    c_dens[i] = (x, y)

lw = 2
fontsize = 14
dperps = [0.1, 0.3, 0.5, 0.7, 0.9]
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, 5))
fig, ax = plt.subplots()

for i in range(2, 7):
    j = i - 2
    label = r"$\mathdefault{D_{\bot}}$ = " + "{:.1f}".format(dperps[j]) + r" $\mathdefault{m^2/s}$"
    ax.plot(c_dens[i][0], c_dens[i][1], color=colors[j], label=label, lw=lw)

ax.set_ylim([0, 4e15])
ax.set_xlabel("Distance from outer target (m)", fontsize=fontsize)
ax.set_ylabel(r"C Density ($\mathdefault{m^{-3}}$)", fontsize=fontsize)
ax.legend(fontsize=fontsize)
ax.grid()

fig.tight_layout()
fig.show()
