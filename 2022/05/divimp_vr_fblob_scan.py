import oedge_plots
import matplotlib.pyplot as plt
from tqdm import tqdm


root = "/Users/zamperini/Documents/d3d_work/divimp_files/blob_test/d3d-167196-blobtest-div-"

# Map names to the values of interest.
vblobs = {"vscan1":200, "vscan2":400, "vscan3":600, "vscan4":800, "vscan5":1000}
fblobs = {"fscan1":0.5e3, "fscan3":1.5e3, "fscan4":2e3, "fscan5":2.5e3, "fscan6":3e3}

# Load objects.
ops = {}
print("Loading runs...")
for i in tqdm(range(1, 7)):
    if i in [1, 2, 3, 4, 5]:
        vpath = root + "vscan{}.nc".format(i)
        ops["vscan{}".format(i)] = oedge_plots.OedgePlots(vpath)
    if i in [1, 3, 4, 5, 6]:
        fpath = root + "fscan{}.nc".format(i)
        ops["fscan{}".format(i)] = oedge_plots.OedgePlots(fpath)

# Get along_ring data.
ring = 25
ar = {}
for k, op in ops.items():
    s, nz = ops[k].along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    ar[k] = {"s":s, "nz":nz}

# Plot comparisions.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), sharey=True)
for k, v in ar.items():
    if k[0] == "v":
        ax1.plot(v["s"], v["nz"], label=vblobs[k])
    else:
        ax2.plot(v["s"], v["nz"], label=fblobs[k])
fig.supxlabel("Distance from inner target (m)")
ax1.set_ylabel("W Density (m-3)")
ax1.legend()
ax2.legend()
ax1.set_title("vblob scan")
ax2.set_title("fblob scan")
ax1.set_yscale("log")
fig.tight_layout()
fig.show()
