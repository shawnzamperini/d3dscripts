# Plot the density profiles as fit from OMFITprofiles, colorcoded by the
# peak frequency as seen on the tile currents.
import netCDF4
import numpy as np
import matplotlib.pyplot as plt


# All LSN, SP on shelf, L-mode. BT = 2T. IP ~ 1.7 MA.
shots = [167196, 170843, 176508, 162766]
freqs = [2030,   2578,   2265,   2317]
bt    = ["r",    "r",    "r",    "f"]

root = "/Users/zamperini/Documents/d3d_work/files/OMFITprofiles_"

data = {}
for shot in shots:
    nc = netCDF4.Dataset("{}{}_FIT.nc".format(root, shot))
    rho = nc["rho"][:]
    ne = nc["n_e"][:].mean(axis=0)
    sep = np.argmin(np.abs(rho - 1))
    data[shot] = {"rho":rho, "ne":ne, "nesep":ne[sep]}

# Color by the frequency.
freqs = np.array(freqs)
norm = plt.Normalize(freqs.min(), freqs.max())
cmap = plt.get_cmap("inferno")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5,4), sharex=True, sharey=True)
ax1.axvline(1.0, color="k")
for i in range(0, len(data)):
    shot = list(data.keys())[i]
    x = data[shot]["rho"]
    y = data[shot]["ne"] / data[shot]["nesep"]
    if bt[i] == "r":
        ax1.plot(x, y, label=shot, color=cmap(norm(freqs[i])))
    elif bt[i] == "f":
        ax2.plot(x, y, label=shot, color=cmap(norm(freqs[i])))
ax1.legend()
ax2.legend()
ax1.set_xlim([1.00, 1.25])
ax1.set_ylim([0.00, 1.2])
ax1.set_xlabel(r"$\rho$", fontsize=16)
ax1.set_ylabel(r"$\mathdefault{n_e/n_{e,sep}}$", fontsize=16)

fig.tight_layout()
fig.show()
