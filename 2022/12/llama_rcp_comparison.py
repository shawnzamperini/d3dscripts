import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import BlobbyFarSOL
from scipy.interpolate import interp1d, splrep, BSpline


path484 = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190484_2.tab"
path485 = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190485_2.tab"
df484 = pd.read_csv(path484, delimiter="\t")
df485 = pd.read_csv(path485, delimiter="\t")

llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190485.npz"
llama = np.load(llama_path)
lrho = llama["rho"]
lpsin = np.sqrt(lrho)  # According to Florian L.
lneut = llama["nDens_LFS"].mean(axis=0)
lneut_err = np.nanmean(llama["nDens_LFS_err"], axis=0)

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
elec = 1.609e-19
md = 2.014 * 1.66e-27  # kg
me = 5.486e-4 * 1.66e-27  # kg
ln_lambda = 15
eps = 8.854e-12  # F/m
all_rcp_rho = np.array([])
all_rcp_nn = np.array([])
all_rcp_midcoll = np.array([])
for plunge in [1, 2]:
    bfs = BlobbyFarSOL.main(190485, plunge)
    min_r = max(bfs.rcp_r.min(), bfs.blob_r.min())
    max_r = min(bfs.rcp_r.max(), bfs.blob_r.max())
    blob_mask = np.logical_and(np.logical_and(bfs.blob_r > min_r, bfs.blob_r < max_r), bfs.neut_dens2 > 0)
    f_rho = interp1d(bfs.rcp_r, bfs.rcp_rho)

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

    all_rcp_rho = np.append(all_rcp_rho, f_rho(bfs.blob_r[blob_mask]))
    all_rcp_nn = np.append(all_rcp_nn, bfs.neut_dens2[blob_mask])
    all_rcp_midcoll = np.append(all_rcp_midcoll, norm_coll)

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
ax.fill_between(lrho, lneut-lneut_err, lneut+lneut_err, color="tab:pink", alpha=0.3)
ax.plot(lrho, lneut, label="LLAMA", color="tab:pink", lw=3)
ax.scatter(all_rcp_rho, all_rcp_nn, color="tab:cyan", s=75, edgecolors="k")
ax.plot(x, yfit, lw=3, color="tab:cyan", label="RCP-Derived")
ax.set_yscale("log")
ax.set_ylim([1e14, 1e18])
ax.set_xlim([0.99, 1.08])
ax.grid(which="both", alpha=0.5)
ax.set_ylabel("Neutral Density (m-3)", fontsize=16)
ax.set_xlabel("Rho", fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=14)
ax.legend(fontsize=14)
fig.tight_layout()
fig.show()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(all_rcp_rho, all_rcp_midcoll)
ax.set_ylabel("Midplane Collisionality", fontsize=16)
ax.set_xlabel("Rho", fontsize=16)
fig.tight_layout()
fig.show()