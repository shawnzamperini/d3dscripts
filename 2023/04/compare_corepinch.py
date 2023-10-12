import oedge_plots
import matplotlib.pyplot as plt


fig, ax = plt.subplots(figsize=(5, 4))
for i in range(0, 6):
    path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-015f-predep-holes-corepinch{}.nc".format(i)
    op = oedge_plots.OedgePlots(path)
    p = op.fake_probe(2.19, 2.36, -0.188, -0.188, "nz", charge="all")
    ax.plot(p["psin"], p["nz"], label=i)

ax.legend()
ax.set_xlabel("Psin")
ax.set_ylabel("nW (m-3)")
ax.set_yscale("log")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.show()