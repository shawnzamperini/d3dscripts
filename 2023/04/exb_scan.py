import oedge_plots
import matplotlib.pyplot as plt


fig, ax = plt.subplots(figsize=(5, 4))
for mult in [0.2, 0.4, 0.6, 0.8, 1.0]:
    path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-013f-{:.1f}exb.nc".format(mult)
    # path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-diff-001-{:.1f}exb.nc".format(mult)
    op = oedge_plots.OedgePlots(path)
    p = op.fake_probe(2.21, 2.36, -0.188, -0.188, "nz", charge="all")
    ax.plot(p["psin"], p["nz"], label=mult)
ax.legend()
ax.set_xlabel("Psin")
ax.set_ylabel("nW (m-3)")
ax.set_yscale("log")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.show()