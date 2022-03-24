# Script to compare the 2D C density results from the current poloidal bumper
# limiters with a continuous toroidal bumper limiter.
from LimPlots import LimPlots
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


root = "/Users/zamperini/Documents/d3d_work/lim_runs/tor_lim_testing/"
angles = [0, 90, 180, 270]
pols = {}
tors = {}
for angle in angles:
    pol_path = "{}167196-pol-bump-tor{}.nc".format(root, angle)
    tor_path = "{}167196-tor-lims-tor{}.nc".format(root, angle)
    pols[angle] = LimPlots(pol_path)
    tors[angle] = LimPlots(tor_path)


pol_nzs = {}
tor_nzs = {}
for angle in angles:
    pol_nzs[angle] = pols[angle].plot_par_rad("nz", 21, showplot=False)
    tor_nzs[angle] = tors[angle].plot_par_rad("nz", 21, showplot=False)

X = pol_nzs[0]["X"]
Y = pol_nzs[0]["Y"]

rmin = -0.055
rmax = -0.05
ymin = 1.0
ymax = 2.0
def plot_mesh(ax, popt, pol, angle):
    if popt == "perc_diff":
        if pol:
            Z = (pol_nzs[angle]["Z"]  - pol_nzs[0]["Z"]) / pol_nzs[0]["Z"]
        else:
            Z = (tor_nzs[angle]["Z"]  - tor_nzs[0]["Z"]) / tor_nzs[0]["Z"]
        vmin = -0.5
        vmax = 0.5
        cmap = "coolwarm"
    elif popt == "nz":
        if pol:
            Z = pol_nzs[angle]["Z"]
        else:
            Z = tor_nzs[angle]["Z"]
        vmin = 5e-5
        vmax = 5e-4
        cmap = "magma"

    ax.pcolormesh(X, Y, Z, shading="auto", cmap=cmap, vmin=vmin,
      vmax=vmax)
    ax.set_xlim([-12, 12])
    if pol:
        bxp = pol_nzs[angle]["pos_bound_x"]
        bxn = pol_nzs[angle]["neg_bound_x"]
        by  = pol_nzs[angle]["bound_y"]
    else:
        bxp = tor_nzs[angle]["pos_bound_x"]
        bxn = tor_nzs[angle]["neg_bound_x"]
        by  = tor_nzs[angle]["bound_y"]
    ax.step(bxp, by, color="r", where="post")
    ax.step(bxn, by, color="r", where="post")
    rect = Rectangle((ymin, rmin), ymax-ymin, rmax-rmin)
    ax.add_patch(rect)

    # Also get the average density in the specified region.
    rmask = np.logical_and(Y>=rmin, Y<=rmax)
    ymask = np.logical_and(X>=ymin, X<=ymax)
    mask = np.logical_and(rmask, ymask)
    avg_nz = Z[mask].mean()
    return {"avg_nz":avg_nz}

fig, axs = plt.subplots(2, 4, figsize=(10, 6))

popt = "nz"
pnzs = []; tnzs = []
pnzs.append(plot_mesh(axs[0][0], popt, True, 0)["avg_nz"])
pnzs.append(plot_mesh(axs[0][1], popt, True, 90)["avg_nz"])
pnzs.append(plot_mesh(axs[1][0], popt, True, 180)["avg_nz"])
pnzs.append(plot_mesh(axs[1][1], popt, True, 270)["avg_nz"])
tnzs.append(plot_mesh(axs[0][2], popt, False, 0)["avg_nz"])
tnzs.append(plot_mesh(axs[0][3], popt, False, 90)["avg_nz"])
tnzs.append(plot_mesh(axs[1][2], popt, False, 180)["avg_nz"])
tnzs.append(plot_mesh(axs[1][3], popt, False, 270)["avg_nz"])

axs[0][0].set_title(r"$\mathdefault{0^{\circ}}$")
axs[0][1].set_title(r"$\mathdefault{90^{\circ}}$")
axs[1][0].set_title(r"$\mathdefault{180^{\circ}}$")
axs[1][1].set_title(r"$\mathdefault{270^{\circ}}$")
axs[0][2].set_title(r"$\mathdefault{0^{\circ}}$")
axs[0][3].set_title(r"$\mathdefault{90^{\circ}}$")
axs[1][2].set_title(r"$\mathdefault{180^{\circ}}$")
axs[1][3].set_title(r"$\mathdefault{270^{\circ}}$")

fig.tight_layout()
fig.show()

pnzs = np.array(pnzs)
tnzs = np.array(tnzs)

# Figure of the avg density at location around the machine.
fig, ax = plt.subplots()
if popt == "perc_diff":
    title = "Percent diff. at R=[{:.3f}, {:.3f}] Y=[{:.1f}, {:.1f}]".format(rmin, rmax, ymin, ymax)
    ax.set_ylabel("Percent difference in nz (%)")
    mult = 100
else:
    title = "Average density at R=[{:.3f}, {:.3f}] Y=[{:.1f}, {:.1f}]".format(rmin, rmax, ymin, ymax)
    ax.set_ylabel("Average impurity density")
    mult = 1
ax.plot(angles, pnzs*mult, "o-", label="Poloidal bumpers", color="tab:red")
ax.plot(angles, tnzs*mult, "o-", label="Toroidal limiters", color="tab:green")
ax.set_xlabel("Toroidal angle")

ax.set_title(title)
ax.grid()
ax.legend()
fig.tight_layout()
fig.show()
