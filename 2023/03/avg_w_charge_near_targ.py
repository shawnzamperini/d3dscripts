import oedge_plots
import matplotlib.pyplot as plt
import numpy as np


blob_nc_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010f.nc"
op = oedge_plots.OedgePlots(blob_nc_path)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)

for charge in range(0, 10):
    d = op.fake_probe(1.42, 1.456, -1.249, -1.249, "nz", charge=charge)
    ax1.plot(d["r"], d["nz"], label=charge)

    if charge == 0:
        wavg1 = np.array(d["nz"]) * charge
        wavg2 = np.array(d["nz"])
    else:
        wavg1 += np.array(d["nz"]) * charge
        wavg2 += np.array(d["nz"])

ax2.plot(d["r"], wavg1 / wavg2)

ax1.axvline(1.404, color="k", linestyle="--")
ax1.axvline(1.454, color="k", linestyle="--")
ax2.axvline(1.404, color="k", linestyle="--")
ax2.axvline(1.454, color="k", linestyle="--")

ax2.axvline(1.423, color="k", linestyle="-")
ax2.axvline(1.423, color="k", linestyle="-")

ax1.legend()
ax1.set_xlabel("R (m)")
ax1.set_ylabel("W Density (m-3)")

ax2.set_ylabel("Average Charge")
ax2.set_xlabel("R (m)")

fig.tight_layout()
fig.show()
