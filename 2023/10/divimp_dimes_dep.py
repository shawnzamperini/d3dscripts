# Script to see what the deposition along a DiMES sample would be from a given DIVIMP simulation.
import oedge_plots
import matplotlib.pyplot as plt


# Load DIVIMP run.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-fluc-002.nc"
op = oedge_plots.OedgePlots(ncpath)

# Pull out deposition along the outer target.
absfac = op.nc["ABSFAC"][:]
deps = op.nc["DEPS"][:].sum(axis=0) * absfac
rp = op.nc["RP"][:]
zp = op.nc["ZP"][:]

# Shelf is at -1.25. May not work for other grids (untested).
shelf_mask = zp == -1.25
rshelf = rp[shelf_mask]
deps_shelf = deps[shelf_mask]

dimes_start = 1.465
dimes_end = 1.515
sp = 1.4225

fig, ax = plt.subplots(figsize=(5, 4))
ax.axvspan(dimes_start, dimes_end, color="grey")
ax.axvline(sp, color="k")
ax.scatter(rshelf, deps_shelf, color="tab:red")
ax.set_xlabel("R (m)", fontsize=14)
ax.set_ylabel("W Areal Deposition (W/m2)", fontsize=14)
ax.set_xlim(1.41, 1.54)
fig.tight_layout()
fig.show()