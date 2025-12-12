# Plot the flows from Flan in cylindrical coordinates to see if they make
# any sense.
import flan_plots
import numpy as np
import matplotlib.pyplot as plt


# Load run
path = "/home/zamp/flandir/sh_outward_v4/sh_outward_v4.nc"
fp = flan_plots.FlanPlots(path)

# Can plot electric field (elec) or ion flow (ion_flow)
dname = "elec"
print("Plotting: {}".format(dname))

# Get X, Y, Z components of flow data
tidx = 4
print("t = {:.2f}".format(fp.nc["time"][tidx]))
uX = fp.nc["{:}_X".format(dname)][tidx]
uY = fp.nc["{:}_Y".format(dname)][tidx]
uZ = fp.nc["{:}_Z".format(dname)][tidx]

# Calculate the radial
uR = np.sqrt(np.square(uX) + np.square(uY))

# Flow parallel to the magnetic field (toroidal, since simple helical) 
# direction
theta = np.arctan2(uY, uX)
uT = -uX * np.sin(theta) + uY * np.cos(theta)

# A radial index in the middle
xidx = int(uX.shape[0] / 2)
print("x = {:.2f}".format(fp.nc["x"][xidx]))

# Index at the midplane (middle y index)
yidx = int(uX.shape[1] / 2)
print("y = {:.2f}".format(fp.nc["y"][yidx]))

# Apply indices, giving the profile along z (the field line)
uR_z = uR[xidx, yidx]
uT_z = uT[xidx, yidx]

# The parallel coordinate at cell center
z = fp.nc["z"][:]

# Plot
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.axhline(0.0, color="k")
ax1.plot(z, uR_z, label="radial")
ax1.plot(z, uT_z, label="parallel")
ax1.legend()
fig.tight_layout()
fig.show()
