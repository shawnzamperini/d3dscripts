import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import BlobbyFarSOL
from scipy.interpolate import interp1d, splrep, BSpline, griddata
from scipy.signal import savgol_filter, medfilt
import pickle
import openadas


shot = 190486

path484 = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190484_2.tab"
path485 = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190485_2.tab"
df484 = pd.read_csv(path484, delimiter="\t")
df485 = pd.read_csv(path485, delimiter="\t")

# Fits for lambda_ne.
mask484 = np.logical_and(df484["R(cm)"] > 224.5, df484["R(cm)"] < 229)
z484 = np.polyfit(df484["R(cm)"][mask484], np.log(df484["Ne(E18 m-3)"][mask484]*1e18), 1)
p484 = np.poly1d(z484)
mask485 = np.logical_and(df485["R(cm)"] > 224.5, df485["R(cm)"] < 232)
z485 = np.polyfit(df485["R(cm)"][mask485], np.log(df485["Ne(E18 m-3)"][mask485]*1e18), 1)
p485 = np.poly1d(z485)
print("190484: {:.2f}".format(-1 / z484[0]))
print("190485: {:.2f}".format(-1 / z485[0]))

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(df484["R(cm)"], df484["Ne(E18 m-3)"]*1e18, label=190484, color="tab:red", lw=3)
ax.plot(df485["R(cm)"], df485["Ne(E18 m-3)"]*1e18, label=190485, color="tab:purple", lw=3)
ax.plot(df484["R(cm)"], np.exp(p484(df484["R(cm)"])), color="tab:red", linestyle="--")
ax.plot(df485["R(cm)"], np.exp(p485(df485["R(cm)"])), color="tab:purple", linestyle="--")
ax.legend()
ax.grid(which="both", alpha=0.5)
ax.set_yscale("log")
ax.set_xlabel("R (cm)", fontsize=14)
ax.set_ylabel(r"$\mathdefault{n_e\ (m^{-3})}$", fontsize=14)
ax.set_ylim([1e18, 1e19])
fig.tight_layout()
fig.show()

boltz = 8.617e-5  # eV/K
elec = 1.6022e-19
md = 2.014 * 1.66e-27  # kg
me = 5.486e-4 * 1.66e-27  # kg
ln_lambda = 15
eps = 8.854e-12  # F/m
data = {}
for shot in [190484, 190485, 190486]:

    # LLAMA data.
    llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_{}_.npz".format(shot)
    llama = np.load(llama_path)
    lpsin = llama["psi_n"]
    lrho = np.square(lpsin)  # According to Florian L.
    lneut = llama["nDens_LFS"]
    lneut_err = llama["nDens_LFS_err"]

    # Get some gfile stuff so we can go from R, Z to psin.
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/{}/{}_3000.pickle".format(shot, shot)
    with open(gfile_path, "rb") as f:
        gfile = pickle.load(f)
    R = gfile["R"]
    Z = gfile["Z"]
    Rs, Zs = np.meshgrid(R, Z)
    psin = gfile["PSIRZ_NORM"]
    Btot = np.sqrt(np.square(gfile["Bt"]) + np.square(gfile["Bz"]))
    all_rcp_rho = np.array([])
    all_rcp_nn = np.array([])
    all_rcp_midcoll = np.array([])
    all_plunge = np.array([])
    all_yinstab = np.array([])
    all_instab_ratio = np.array([])
    all_norm_vr = np.array([])
    for plunge in [1, 2]:
        bfs = BlobbyFarSOL.main(shot, plunge, showplot=False, temod=0.50)
        min_r = max(bfs.rcp_r.min(), bfs.blob_r.min())
        max_r = min(bfs.rcp_r.max(), bfs.blob_r.max())
        blob_mask = np.logical_and(np.logical_and(bfs.blob_r > min_r, bfs.blob_r < max_r), bfs.neut_dens2 > 0)
        f_rho = interp1d(bfs.rcp_r, bfs.rcp_rho)

        rcp_coord = zip(bfs.blob_r / 100, np.full(len(bfs.blob_r), -0.185))
        rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))

        # Midplane collisionality.
        # f_conn = interp1d(bfs.conns_r, bfs.conns_l_otf + bfs.conns_l_itf)
        # nu_ei_hat = 1.33e5 * bfs.blob_ne[blob_mask] * 1e-20 / np.power(bfs.blob_te[blob_mask] * 1e-3, 3 / 2)
        # nu_ei = 1.4 * (5.486e-4 / 2.014) * nu_ei_hat  # 1.4 * (me/mi) * nu_hat
        # elec_gyrofreq = 1.76e11 * bfs.blob_b[blob_mask]  # Eq. 8.20 in Friedberg for the electron gyrofreq.
        # ion_gyrorad = 6.46e-3 * np.sqrt(bfs.blob_te[blob_mask] * 1e-3) / bfs.blob_b[blob_mask]
        # norm_coll = nu_ei * f_conn(bfs.blob_r[blob_mask]) / (elec_gyrofreq * ion_gyrorad)
        # norm_coll2 = 2.75e-14 * f_conn(bfs.blob_r[blob_mask]) * bfs.blob_ne[blob_mask] / np.square(bfs.blob_te[blob_mask])

        # Need to calculate the normalized collisionality and blob size. First calculate collisionality.
        f_conn_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
        vi = np.sqrt(3 * elec * bfs.blob_te[blob_mask] / md)   # m/s
        vte = np.sqrt(2 * elec * bfs.blob_te[blob_mask] / md)  # m/s
        nu_ei = (1 / (4 * np.pi) * bfs.blob_ne[blob_mask] * elec**4 / (eps**2 * me * md) * ln_lambda) * \
                (1 / (np.power(vi, 3) + 1.3 * np.power(vte, 3)))   # 1/s
        elec_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / me
        ion_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / md
        ion_gyrorad = np.sqrt(2 * 2.014 * 1.66e-27 * elec * bfs.blob_te[blob_mask]) / (elec * bfs.blob_b[blob_mask])
        norm_coll = nu_ei * f_conn_otf(bfs.blob_r[blob_mask]) / (elec_gyrofreq * ion_gyrorad)

        # Now calculate the normalized radius.
        ab = bfs.blob_drad[blob_mask] / 2 / 100  # cm to m
        R = bfs.blob_r[blob_mask] / 100  # cm to m
        L = f_conn_otf(bfs.blob_r[blob_mask])
        ahat = ab * np.power(R, 1 / 5) / (np.power(L, 2 / 5) * np.power(ion_gyrorad, 4 / 5))
        norm_rad = np.power(ahat, 5 / 2)

        # Now the normalized velocity.
        vr = bfs.blob_vr[blob_mask]
        cs = bfs.blob_cs[blob_mask]
        norm_vr = vr / (np.power(2 * L * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs)

        # Calculate the blob density above the background value.
        f_rcp_ne = interp1d(bfs.rcp_r, bfs.rcp_ne)
        deltn_n = (bfs.blob_ne[blob_mask] - f_rcp_ne(bfs.blob_r[blob_mask])) / f_rcp_ne(bfs.blob_r[blob_mask])
        print("deltn_n")
        print(deltn_n)
        print()

        # Calculate the ion-neutral collision (ionization) rate.
        nu_iz = bfs.neut_dens2[blob_mask] * bfs.blob_ioniz_rate_coeffs[blob_mask]

        # Now calculate the secondary instability growth rate, Eq. 5.1.13 from Theiler's thesis.
        blob_rzs = zip(R, np.full(len(R), -0.185))
        B = griddata((Rs.flatten(), Zs.flatten()), Btot.flatten(), list(blob_rzs))
        te = bfs.blob_te[blob_mask]
        larmor = np.sqrt(te * md * elec) / (elec * B)
        yint = np.sqrt(2 / (R * ab)) * cs
        instab_ratio = np.sqrt(2 * ab / R) * cs * deltn_n / vr - (1 / (np.square(larmor) * L)) * np.sqrt(R / 2) \
                       * np.power(ab, 5/2) - nu_iz * np.sqrt(R * ab) / (np.sqrt(2) * cs)
        yinstab = instab_ratio * yint
        print("instab_ratio")
        print(instab_ratio)
        print()

        all_rcp_rho = np.append(all_rcp_rho, rcp_psin[blob_mask])
        all_rcp_nn = np.append(all_rcp_nn, bfs.neut_dens2[blob_mask])
        all_rcp_midcoll = np.append(all_rcp_midcoll, norm_coll)
        all_plunge = np.append(all_plunge, np.full(len(rcp_psin[blob_mask]), plunge))
        all_yinstab = np.append(all_yinstab, yinstab)
        all_instab_ratio = np.append(all_instab_ratio, instab_ratio)
        all_norm_vr = np.append(all_norm_vr, norm_vr)
    sort_idx = np.argsort(all_rcp_rho)
    all_rcp_rho = all_rcp_rho[sort_idx]
    all_rcp_nn = all_rcp_nn[sort_idx]
    all_rcp_midcoll = all_rcp_midcoll[sort_idx]
    all_plunge = all_plunge[sort_idx]
    all_yinstab = all_yinstab[sort_idx]
    all_instab_ratio = all_instab_ratio[sort_idx]
    all_norm_vr = all_norm_vr[sort_idx]
    data[shot] = {"psin": all_rcp_rho, "nn": all_rcp_nn, "plunge":all_plunge, "lpsin": lpsin, "lneut": lneut,
                  "yinstab":all_yinstab, "instab_ratio":all_instab_ratio, "norm_vr":all_norm_vr}

# Spline fit.
sort_idx = np.argsort(all_rcp_rho)
x = all_rcp_rho[sort_idx][2:]
y = all_rcp_nn[sort_idx][2:] / 1e17
knot_numbers = 2
x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
q_knots = np.quantile(x, x_new)
t,c,k = splrep(x, y, s=50)
yfit = BSpline(t, c, k)(x) * 1e17


fig, ax = plt.subplots(figsize=(6, 5))
ax.axvline(1.0, color="k", linestyle="--")
# ax.fill_between(lpsin, lneut-lneut_err, lneut+lneut_err, color="tab:pink", alpha=0.3)
ax.plot(lpsin, lneut, label="LLAMA", color="tab:pink", lw=3, zorder=15)
ax.scatter(all_rcp_rho, all_rcp_nn, color="tab:cyan", s=75, edgecolors="k", zorder=25)
ax.plot(x, yfit, lw=3, color="tab:cyan", label="RCP-Derived", zorder=35)
ax.set_yscale("log")
ax.set_ylim([1e14, 1e18])
ax.set_xlim([0.98, 1.20])
ax.grid(which="both", alpha=0.5)
ax.set_ylabel("Neutral Density (m-3)", fontsize=16)
ax.set_xlabel("Psin", fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=14)
ax.legend(fontsize=14)
fig.tight_layout()
fig.show()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(all_rcp_rho, all_rcp_midcoll)
ax.set_ylabel("Midplane Collisionality", fontsize=16)
ax.set_xlabel("Psin", fontsize=16)
fig.tight_layout()
fig.show()

# Load a DIVIMP background from a SAS-VW discharge to get a sense of the order of magnitude difference in neutral
# densities between the two locations.
import oedge_plots
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190484/d3d-190484-bkg-005-outgas-smooth.nc"
op = oedge_plots.OedgePlots(ncpath)
div_llama = op.fake_probe(1.92, 2.05, -0.77, -0.77, data="neut_dens", plot="psin", show_plot=False)
div_rcp = op.fake_probe(2.23, 2.33, -0.185, -0.185, data="neut_dens", plot="psin", show_plot=False)

fig, ax1 = plt.subplots(figsize=(5, 4))
for shot in [190484, 190485, 190486]:
    if shot == 190484:
        color = "tab:purple"
    elif shot == 190485:
        color = "tab:green"
    elif shot == 190486:
        color = "tab:red"

    window = 51
    lneutfit = savgol_filter(medfilt(data[shot]["lneut"], window), window, 3)
    nnfit = savgol_filter(medfilt(data[shot]["nn"], 5), 5, 3)
    ax1.axvline(1.0, linestyle="--", color="k")
    ax1.plot(data[shot]["lpsin"], lneutfit, color=color, lw=3)
    # ax1.scatter(data[shot]["psin"], data[shot]["nn"], color=color)
    ax1.plot(data[shot]["psin"], nnfit, color="k", lw=4)
    ax1.plot(data[shot]["psin"], nnfit, color=color, lw=3)
    ax1.set_xlabel("Psin", fontsize=14)
    ax1.set_ylabel("Neutral Density", fontsize=14)
    ax1.set_yscale("log")
    ax1.set_xlim([0.98, 1.15])
    ax1.set_ylim([1e14, 1e18])
    ax1.grid(which="both", alpha=0.3)
custom_lines = [Line2D([0], [0], color="tab:purple", lw=3),
                Line2D([0], [0], color="tab:green", lw=3),
                Line2D([0], [0], color="tab:red", lw=3)]
ax1.legend(custom_lines, ["190484", "190485", "190486"], fontsize=14, framealpha=1.0)
fig.tight_layout()
fig.show()

window = 51
lneutfit = savgol_filter(medfilt(data[190484]["lneut"], window), window, 3)
nnfit = savgol_filter(medfilt(data[190484]["nn"], 5), 5, 3)

fig, ax2 = plt.subplots(figsize=(5,4))
ax2.axvline(1.0, color="k", linestyle="--")
ax2.plot(data[190484]["lpsin"], lneutfit, color="tab:red", lw=3)
ax2.scatter(data[190484]["psin"], data[190484]["nn"], color="tab:purple", zorder=35, edgecolors="k", s=150, marker="*")
ax2.plot(div_llama["psin"], div_llama["neut_dens"], lw=5, color="k")
ax2.plot(div_llama["psin"], div_llama["neut_dens"], label="LLAMA", lw=3, color="tab:red")
ax2.plot(div_rcp["psin"], div_rcp["neut_dens"], lw=5, color="k")
ax2.plot(div_rcp["psin"], div_rcp["neut_dens"], label="RCP", lw=3, color="tab:purple")
ax2.legend()
ax2.set_yscale("log")
ax2.grid(alpha=0.3, which="both")
ax2.set_xlabel("Psin", fontsize=14)
ax2.set_ylabel("Neutral Density (m-3)", fontsize=14)
ax2.set_title("#190484", fontsize=14)
ax2.set_xlim([0.97, 1.15])
ax2.set_ylim([5e14, 1e17])
fig.tight_layout()
fig.show()

# Plot involving the secondary linear growth rates.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
ax1.axvline(1.0, color="k", linestyle="--")
ax1.axhline(0, color="k", linestyle="-")
for shot in [190484, 190485, 190486]:
    if shot == 190484:
        color = "tab:purple"
        window = 5
    elif shot == 190485:
        color = "tab:green"
        window = 3
    elif shot == 190486:
        color = "tab:red"
        window = 7
    ax1.plot(data[shot]["psin"], data[shot]["instab_ratio"], color=color, lw=3, label=shot)

    ax2.scatter(data[shot]["instab_ratio"], data[shot]["norm_vr"], label=shot, zorder=15, color=color, edgecolors="k")

#ax1.set_yscale("log")
ax1.set_xlabel("Psin", fontsize=14)
ax1.set_ylabel(r"$\gamma_{inst}$ (s)", fontsize=14)
ax1.grid(which="both", alpha=0.3)
ax1.set_xlim([0.97, 1.15])
ax1.legend(loc="lower left")
ax2.grid(zorder=5, alpha=0.3)
fig.tight_layout()
fig.show()