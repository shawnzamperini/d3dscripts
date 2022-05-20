# Plot output from slab_blob.
import matplotlib.pyplot as plt
import pickle


with open("slab_blob_output.pickle", "rb") as f:
    results = pickle.load(f)


X = results["X"]
Y = results["Y"]
nb = results["blob_dens"]
targ1 = results["targ1"]
targ2 = results["targ2"]
targ1_dep = results["targ1_dep"]
targ2_dep = results["targ2_dep"]
te = results["te"]
fblob_sep = results["fblob_sep"]
nblobs = results["nblobs"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,8))

ax1.pcolormesh(X, Y, nb, shading="auto")
ax1.set_xlabel("Radial (m)")
ax1.set_ylabel("Parallel (m)")

ax2.plot(targ1, targ1_dep/nblobs*fblob_sep)

fig.tight_layout()
fig.show()
