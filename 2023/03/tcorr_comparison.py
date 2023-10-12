# Script to make a quick plot comparing a scan in the correlation time.
import oedge_plots
import matplotlib.pyplot as plt


tcorrs = [100, 50, 10, 200, 500, 5, 20]  # us
keys = ["", "a", "b", "c", "d", "f", "g"]
ops = []
nzs = []
for i in range(0, len(keys)):
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010{}.nc".format(keys[i])
    op = oedge_plots.OedgePlots(ncpath)
    ops.append(op)

    s, nz = op.along_ring(91, "DDLIMS", charge="all", plot_it=False)
    nzs.append(nz)

# The case from my FST paper that had a decent match in the 3DLIM sim with it.
fst_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mrc-shifted-nodrift-2-copy.nc"
fst = oedge_plots.OedgePlots(fst_path)
s_fst, nz_fst = fst.along_ring(91, "DDLIMS", charge="all", plot_it=False)

fig, ax = plt.subplots(figsize=(5, 4))
for i in range(0, len(keys)):
    ax.plot(s, nzs[i], label=tcorrs[i])
ax.plot(s_fst, nz_fst, color="r", lw=3, label="FST")
ax.legend()
ax.set_xlabel("S (m)")
ax.set_ylabel("W Density (m-3)")
fig.tight_layout()
fig.show()
