# Make a plot for my DIVIMP blobby paper summarizing the RCP measurements for a given discharge.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import pickle


# Load data.
shot = 167195; plunge = 2
shot = 190411; plunge = 2
if shot == 167196:
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
    mafot_m1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp4500_conn_-1.dat"
    mafot_p1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp4500_conn_+1.dat"
elif shot == 190411:
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190411/190411_3000.pickle"  # Same shape for the most part.
    mafot_m1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190411/lam_rcp3000_conn_-1.dat"
    mafot_p1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190411/lam_rcp3000_conn_+1.dat"

rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP{}_{}.tab".format(shot, plunge)
rcp = pd.read_csv(rcp_path, delimiter="\t")
blob_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_{}_{}.tab".format(shot, plunge)
blob = pd.read_csv(blob_path, delimiter="\t")
timestep = (blob["Time(ms)"].iloc[1] - blob["Time(ms)"].iloc[0]) / 1000  # ms to s
fblob = blob["Npeaks"].values / timestep

# Get data on psin.
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
rcp_coord = zip(rcp["R(cm)"] / 100, np.full(len(rcp), -0.185))
blob_coord = zip(blob["R(cm)"] / 100, np.full(len(blob), -0.185))
rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))
blob_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(blob_coord))

if shot == 167195:
    rcp_mask = rcp_psin > 1.077
    blob_mask = blob_psin > 0
elif shot == 190411:
    rcp_mask = rcp_psin < 1.21
    blob_mask = blob_psin > 1.05

# Load some connection length data.
def load_mafot(itf_path, otf_path):
    print("Loading MAFOT data...")
    columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
               "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
    mafot_itf = pd.read_csv(itf_path, skiprows=52, names=columns, delimiter="\t")
    mafot_otf = pd.read_csv(otf_path, skiprows=52, names=columns, delimiter="\t")
    conns_r = mafot_itf["R (m)"]
    conns_l_itf = mafot_itf["Lconn (km)"].values * 1000  # km to m
    conns_l_otf = mafot_otf["Lconn (km)"].values * 1000
    return {"r":conns_r, "litf":conns_l_itf, "lotf":conns_l_otf}


mafot = load_mafot(mafot_m1_path, mafot_p1_path)
mafot_coord = zip(mafot["r"], np.full(len(mafot["r"]), -0.185))
mafot_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(mafot_coord))

lw = 2
color = "tab:red"
fig, axs = plt.subplots(4, 2, sharex=True, figsize=(5, 6), zorder=5)
axs = np.array(axs).flatten()
# fig.text(0.5, 0.00, r"$\mathdefault{\psi_n}$", ha="center")
# fig = plt.figure(figsize=(5, 6))
# big_ax = fig.add_subplot(111, frameon=False)
# big_ax.spines['top'].set_color('none')
# big_ax.spines['bottom'].set_color('none')
# big_ax.spines['left'].set_color('none')
# big_ax.spines['right'].set_color('none')
# big_ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
# big_ax.set_xlabel(r"$\mathdefault{\psi_n}$")
# axs = []
# for i in range(1, 9):
#     axs.append(fig.add_subplot(int("42{}".format(i))))

axs[0].plot(rcp_psin[rcp_mask], rcp["Te(eV)"][rcp_mask], lw=lw, color=color, label=r"$\mathdefault{T_e}$ (eV)")

axs[2].plot(rcp_psin[rcp_mask], rcp["Ne(E18 m-3)"][rcp_mask]*1e18, lw=lw, color=color,
            label=r"$\mathdefault{n_e\ (m^{-3})}$")
if shot == 167195:
    axs[2].set_ylim([1e18, 2e19])
elif shot == 190411:
    axs[2].set_ylim([8e17, 2e19])
axs[2].set_yscale("log")

axs[4].plot(rcp_psin[rcp_mask], rcp["Q_par(MW/m^2)"][rcp_mask], lw=lw, color=color,
            label=r"$\mathdefault{q_{\|\|}\ (MW/m^2)}$")
if shot == 167195:
    axs[4].set_ylim([0.01, 10])
elif shot == 190411:
    axs[4].set_ylim([1e-2, 2e1])
axs[4].set_yscale("log")

axs[6].axhline(0, color="k")
axs[6].plot(rcp_psin[rcp_mask], rcp["Machn"][rcp_mask], lw=lw, color=color, label="Mach")

axs[1].plot(blob_psin[blob_mask], fblob[blob_mask], lw=lw, color=color, label=r"$\mathdefault{f_{blob}\ (Hz)}$")

# axs[3].plot(blob_psin, blob["Epol(V/m)"], lw=lw, color=color, label=r"$\mathdefault{E_{\theta}\ (V/m)}$")

axs[3].plot(blob_psin[blob_mask], blob["Vr(m/s)"][blob_mask], lw=lw, color=color, label=r"$\mathdefault{v_r}$ (m/s)")

axs[5].plot(blob_psin[blob_mask], blob["D_rad(cm)"][blob_mask], lw=lw, color=color, label=r"$\mathdefault{D_{rad}}$ (cm)")

axs[7].plot(mafot_psin, mafot["litf"], lw=lw, color="tab:purple", label=r"$\mathdefault{L_{ITF}}$")
axs[7].plot(mafot_psin, mafot["lotf"], lw=lw, color="tab:purple", linestyle="--", label=r"$\mathdefault{L_{OTF}}$")
axs[7].set_yscale("log")
axs[7].set_ylim([0.5, 100])

axs[6].set_xlabel(r"$\mathdefault{\psi_N}$")
axs[7].set_xlabel(r"$\mathdefault{\psi_N}$")
axs[0].set_title("Standard Analysis")
axs[1].set_title("Blob Analysis")

i = 0
labels = ["a", "e", "b", "f", "c", "g", "d", "h"]
for ax in axs:
    if ax in [axs[2]]:
        ax.grid(alpha=0.3, which="both")
    else:
        ax.grid(alpha=0.3)

    if shot == 167195:
        if ax in [axs[2], axs[4], axs[6], axs[7]]:
            ax.legend(loc="lower left", ncol=2)
        else:
            ax.legend(loc="upper right")
    elif shot == 190411:
        if ax in [axs[1], axs[3], axs[5], axs[4], axs[6], axs[7]]:
            ax.legend(loc="lower left", ncol=2)
        else:
            ax.legend(loc="upper right")


    if shot == 167196:
        ax.set_xlim([1.05, 1.3])
    elif shot == 190411:
        ax.set_xlim([0.99, 1.22])

    ax.text(0.03, 0.85, labels[i]+")", transform=ax.transAxes)
    i += 1
fig.tight_layout()
fig.show()