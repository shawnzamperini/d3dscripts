import postgkyl
import matplotlib.pyplot as plt
import numpy as np


# Result root
#root = "/home/zamp/gkyldir/NT-current/gk_d3d_negD_iwl_3x2v"
root = "/home/zamp/gkyldir/NT-aug25-60x60x16res/gk_d3d_iwl_adapt_source_3x2v_p1"

# Load covariant components of magnetic field
gdata_b_i = postgkyl.data.GData(root + "-b_i.gkyl")
interp_b_i = postgkyl.data.GInterpModal(gdata_b_i, 1, "ms")
b_i = []
for i in range(3):
	grid, tmp = interp_b_i.interpolate(comp=i)
	b_i.append(tmp[:,:,:,0])

# Reciprocal basis vectors
gdata_dzdx = postgkyl.data.GData(root + "-dzdx.gkyl")
interp_dzdx = postgkyl.data.GInterpModal(gdata_dzdx, 1, "ms")
dzdx = []
for i in range(9):
	grid, tmp = interp_dzdx.interpolate(i)
	dzdx.append(tmp[:,:,:,0])

# Cartesian components of magnetic field (unit?) vector
# Way it's coded into Flan
bX = b_i[0] * dzdx[0] + b_i[1] * dzdx[3] + b_i[2] * dzdx[6]
bY = b_i[0] * dzdx[1] + b_i[1] * dzdx[4] + b_i[2] * dzdx[7]
bZ = b_i[0] * dzdx[2] + b_i[1] * dzdx[5] + b_i[2] * dzdx[8]

# Alternative
#bX = b_i[0] * dzdx[0] + b_i[1] * dzdx[1] + b_i[2] * dzdx[2]
#bY = b_i[0] * dzdx[3] + b_i[1] * dzdx[4] + b_i[2] * dzdx[5]
#bZ = b_i[0] * dzdx[6] + b_i[1] * dzdx[7] + b_i[2] * dzdx[8]

# Prepare for plotting
yidx = 32 # Midpoint of y range
x = grid[0]
z = grid[2]
bX_slice = bX[:, yidx, :].T
bY_slice = bY[:, yidx, :].T
bZ_slice = bZ[:, yidx, :].T
X, Z = np.meshgrid(x, z)
minval = np.min([bX_slice.min(), bY_slice.min(), bZ_slice.min()])
maxval = np.max([bX_slice.max(), bY_slice.max(), bZ_slice.max()])
clim = np.abs([minval, maxval]).max()

# Plot
fig, axes = plt.subplots(1, 3, figsize=(12,4), constrained_layout=True)
fig.suptitle("yidx = {}".format(yidx))
for ax, bcomp, comp in zip(axes, [bX_slice, bY_slice, bZ_slice], ["BX", "BY", "BZ"]):
	mesh = ax.pcolormesh(X, Z, bcomp, cmap="coolwarm", vmin=-clim, vmax=clim)
	ax.set_xlabel("x")
	ax.set_ylabel("z")
	ax.set_title(comp)
cbar = fig.colorbar(mesh, ax=axes, orientation='vertical', shrink=0.8)
cbar.set_label("B (T)")
fig.show()
