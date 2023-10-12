import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import griddata
import oedge_plots

sys.path.append("../../2022/12")
import BlobbyFarSOL

data = {}
for shot in [190484, 190485, 190486]:
    print(shot)

    if shot in [190442]:
        llama_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/LLAMA_190422.npz"
    if shot == 190484:
        llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190484_.npz"
        gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle"
    if shot == 190485:
        llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190485_.npz"
        gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190485/190485_3000.pickle"
    if shot == 190486:
        llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190486_.npz"
        gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190486/190486_3000.pickle"

    llama = np.load(llama_path)
    lrho = np.square(llama["psi_n"])
    lneut = llama["nDens_LFS"]
    lneut_err = llama["nDens_LFS_err"]
    lion = llama["ion_LFS"]
    lion_err = llama["ion_LFS_err"]
    lpsin = llama["psi_n"]

    with open(gfile_path, "rb") as f:
        gfile = pickle.load(f)

    # Get some gfile stuff so we can go from R, Z to psin.
    R = gfile["R"]
    Z = gfile["Z"]
    Rs, Zs = np.meshgrid(R, Z)
    psin = gfile["PSIRZ_NORM"]

    # Load the bfs object as needed.
    bfs_data = {}
    for plunge in [1, 2]:
        bfs = BlobbyFarSOL.main(shot, plunge, showplot=False, temod=1.0)

        rcp_coord = zip(bfs.rcp_r / 100, np.full(len(bfs.rcp_r), -0.185))
        rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))

        bfs_data[plunge] = {"r": bfs.rcp_r, "psin": rcp_psin, "nn": bfs.neut_dens7, "ne":bfs.rcp_ne, 
                            "blob_gamma":bfs.blob_rad_flux, "blob_r":bfs.blob_r}

    data[shot] = {"rcp_psin1": bfs_data[1]["psin"], "rcp_psin2": bfs_data[2]["psin"], "rcp_nn1": bfs_data[1]["nn"],
                  "rcp_nn2": bfs_data[2]["nn"], "llama_psin": lpsin, "llama_neut": lneut, "llama_neut_err": lneut_err,
                  "rcp_ne1":bfs_data[1]["ne"], "rcp_ne2":bfs_data[2]["ne"], "rcp_blob_flux1":bfs_data[1]["blob_gamma"],
                  "rcp_blob_flux2":bfs_data[2]["blob_gamma"], "blob_r1":bfs_data[1]["blob_r"], "blob_r2":bfs_data[2]["blob_r"]}

err_thresh = 0.15
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 4), sharex=True, sharey=True)
for shot, ax in ((190484, ax1), (190485, ax2), (190486, ax3)):

    # For LLAMA data, remove data above a set error threshold.
    lx = data[shot]["llama_psin"]
    ly = data[shot]["llama_neut"]
    lyerr = data[shot]["llama_neut_err"]
    mask = lyerr / ly < err_thresh
    lx = lx[mask]
    ly = ly[mask]

    ax.axvline(1.0, color="k", linestyle="--", zorder=7)
    ax.scatter(lx, ly, s=25, edgecolors="k", color="tab:cyan", label="LLAMA", zorder=10)
    # ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], lw=3, color="k")
    # ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], lw=2, color="tab:red", label="Simple SOL")
    # ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_nn2"], lw=3, color="k")
    # ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_nn2"], lw=2, color="tab:red")
    ax.scatter(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], color="tab:red", label="TSS (Plunge #1)", edgecolors="k",
               marker="*", s=85, zorder=15)
    ax.scatter(data[shot]["rcp_psin2"], data[shot]["rcp_nn2"], color="tab:purple", edgecolors="k", marker="*", s=85,
               zorder=15, label="TSS (Plunge #2)")
    ax.grid(alpha=0.3, zorder=5)
    ax.set_xlabel(r"$\mathdefault{\psi_N}$", fontsize=12)
    ax.set_title("#{}".format(shot), fontsize=12)

    # Density on the secondary axis.
    #ax22 = ax.twinx()
    # ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_ne1"], color="k", zorder=12, lw=3)
    # ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_ne2"], color="k", zorder=12, lw=3)
    # ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_ne1"], color="tab:red", zorder=12, lw=2)
    # ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_ne2"], color="tab:purple", zorder=12, lw=2)
    # ax.set_yscale("log")
    # ax.set_ylim([7e17, 1e19])

ax1.set_yscale("log")
ax1.legend()
ax1.set_ylabel(r"Neutral Density ($\mathdefault{m^{-3}}$)", fontsize=12)
ax1.set_xlim([0.95, 1.25])
ax1.set_ylim([1e14, 1e21])
#ax22.set_ylabel(r"$\mathdefault{n_e\ (m^{-3})}$", fontsize=12)
fig.tight_layout()
fig.show()


# Additional OEDGE run where we've included wall outgassing to simulate the neutrals.
# Load background.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190484/d3d-190484-bkg-005-outgas-smooth.nc"
op = oedge_plots.OedgePlots(ncpath)

# Get neutral densities at each location.
div_llama = op.fake_probe(1.93, 2.05, -0.77, -0.77, data="neut_dens", plot="psin", show_plot=False)
div_rcp   = op.fake_probe(2.23, 2.33, -0.185, -0.185, data="neut_dens", plot="psin", show_plot=False)

# Create plot.
fig, ax = plt.subplots(figsize=(4, 3))

# For LLAMA data, remove data above a set error threshold.
for shot in [190484, 190485, 190486]:
    lx = data[shot]["llama_psin"]
    ly = data[shot]["llama_neut"]
    lyerr = data[shot]["llama_neut_err"]
    mask = lyerr / ly < err_thresh
    lx = lx[mask]
    ly = ly[mask]
    if shot == 190484:
        ax.scatter(lx, ly, s=25, edgecolors="k", color="tab:cyan", label="LLAMA")
    else:
        ax.scatter(lx, ly, s=25, edgecolors="k", color="tab:cyan")

ax.axvline(1.0, color="k", linestyle="--")

for shot in [190484, 190485, 190486]:
    # ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], lw=3, color="k")
    # ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], lw=2, color="tab:red", label="Simple SOL")
    # ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_nn2"], lw=3, color="k")
    # ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_nn2"], lw=2, color="tab:red")
    if shot == 190484:
        ax.scatter(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], color="tab:red", label="Simple SOL", edgecolors="k",
                   marker="*", s=75, zorder=15)
    else:
        ax.scatter(data[shot]["rcp_psin1"], data[shot]["rcp_nn1"], color="tab:red", edgecolors="k",
                   marker="*", s=75, zorder=15)
    ax.scatter(data[shot]["rcp_psin2"], data[shot]["rcp_nn2"], color="tab:red", edgecolors="k", marker="*", s=75,
               zorder=15)
#ax.plot(div_llama["psin"], div_llama["neut_dens"], color="tab:cyan")
#ax.plot(div_rcp["psin"], div_rcp["neut_dens"], color="tab:red")

ax.set_yscale("log")
ax.legend()
ax.grid(alpha=0.3, which="both", zorder=5)
ax.set_xlabel(r"$\mathdefault{\psi_N}$")
ax.set_ylabel(r"Neutral Density ($\mathdefault{m^{-3}}$)", fontsize=12)
ax.set_xlim([0.95, 1.25])
ax.set_ylim([1e14, 1e21])
fig.tight_layout()
fig.show()

shot = 190485
fig, ax = plt.subplots(figsize=(5, 4))
ax.axvline(1.0, color="k", linestyle="--")
ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_ne1"], color="k", zorder=12, lw=3)
ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_ne2"], color="k", zorder=12, lw=3)
ax.plot(data[shot]["rcp_psin1"], data[shot]["rcp_ne1"], color="tab:red", zorder=12, lw=2, label="Plunge #1")
ax.plot(data[shot]["rcp_psin2"], data[shot]["rcp_ne2"], color="tab:purple", zorder=12, lw=2, label="Plunge #2")
ax.set_yscale("log")
ax.legend()
ax.grid(alpha=0.3, zorder=5, which="both")
ax.set_xlabel(r"$\mathdefault{\psi_N}$", fontsize=12)
ax.set_ylabel(r"$\mathdefault{n_e\ (m^{-3})}$", fontsize=12)
ax.set_ylim([1e18, 1e19])
fig.tight_layout()
fig.show()

x1 = data[shot]["blob_r1"]
x2 = data[shot]["blob_r2"]
lny1 = np.log(data[shot]["rcp_blob_flux1"])
lny2 = np.log(data[shot]["rcp_blob_flux2"])
z1 = np.polyfit(x1, lny1, 1)
z2 = np.polyfit(x2, lny2, 1)
fity1 = np.exp(z1[0] * x1 + z1[1])
fity2 = np.exp(z2[0] * x2 + z2[1])

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(data[shot]["blob_r1"], data[shot]["rcp_blob_flux1"], color="tab:red", zorder=12, edgecolors="k", alpha=0.3, label="Plunge #1")
ax.scatter(data[shot]["blob_r2"], data[shot]["rcp_blob_flux2"], color="tab:purple", zorder=12, edgecolors="k", alpha=0.3, label="Plunge #2")
ax.plot(x1, fity1, color="k", lw=3, zorder=15)
ax.plot(x2, fity2, color="k", lw=3, zorder=15)
ax.plot(x1, fity1, color="tab:red", lw=2, zorder=15)
ax.plot(x2, fity2, color="tab:purple", lw=2, zorder=15)
ax.set_yscale("log")
ax.legend()
ax.grid(alpha=0.3, zorder=5, which="both")
ax.set_xlabel("R (cm)", fontsize=12)
ax.set_ylabel(r"$\mathdefault{\tilde{\bar{\Gamma_r}}\ (m^{-2}s^{-1})}$", fontsize=12)
#ax.set_ylim([1e18, 1e19])
fig.tight_layout()
fig.show()