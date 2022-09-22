# This is a weird script to plot a field line trace for 167196 with the wall,
# with the goal of seeing if an RCP trace could hit the low power helicon.
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# First load the EFIT wall and stuff.
gfile_path = "/Users/zamperini/Documents/d3d_work/lwt/toroidal_limiter/167196_3500.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
vesr = gfile["RLIM"]
vesz = gfile["ZLIM"]
sepr = gfile["RBBBS"]
sepz = gfile["ZBBBS"]
gR, gZ = np.meshgrid(gfile["R"], gfile["Z"])
psin = gfile["PSIRZ_NORM"]

# Now load the 3D field line trace from a specified location.
m1path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/struct_lph_m1.dat"
p1path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/struct_lph_p1.dat"
m1 = pd.read_csv(m1path, skiprows=1, names=["X", "Y", "Z", "R", "phi", "psi", "L", "dpsidLc"], delimiter="\t")
p1 = pd.read_csv(p1path, skiprows=1, names=["X", "Y", "Z", "R", "phi", "psi", "L", "dpsidLc"], delimiter="\t")
m1["degree"] = np.degrees(m1["phi"])
p1["degree"] = np.degrees(p1["phi"])

# Pull out the R, Z of where on the vessel the field line terminates for each direction.
rend1 = m1["R"].values[0]
zend1 = m1["Z"].values[0]
rend2 = p1["R"].values[-1]
zend2 = p1["Z"].values[-1]
rstart = 2.31
zstart = -0.185
psi = p1["psi"].mean()

fig, ax1 = plt.subplots()
ax1.set_aspect("equal")
ax1.plot(vesr, vesz, color="k")
ax1.plot(sepr, sepz, color="k")
ax1.contour(gR, gZ, psin, levels=[psi])
ax1.scatter(rstart, zstart, marker="x", color="k", s=50)
ax1.scatter(rend1, zend1, marker="x", color="r", s=50)
ax1.scatter(rend2, zend2, marker="x", color="r", s=50)
fig.tight_layout()
fig.show()
