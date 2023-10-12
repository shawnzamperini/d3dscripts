# Using Tess's full equations, let's make some of the plots we discussed.
import sys
sys.path.append("../../2022/12")
import BlobbyFarSOL
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d, griddata
import pickle


# Inputs to model.
ti_mult = 1.0
shot_plot = 167195
plunge_plot = 2

# Constants.
boltz = 8.617e-5  # eV/K
elec = 1.602e-19
md = 2.014 * 1.66e-27  # kg
md_ev = 2.014 * 931.49e6 / np.square(3e8)  # eV s2 / m2
me = 5.486e-4 * 1.66e-27  # kg
ln_lambda = 15
eps = 8.854e-12  # F/m

vrs = np.array([]); vrs_pred = np.array([]); abs = np.array([]); norm_vrs = np.array([]); norm_vrs_pred = np.array([])
ahats = np.array([]); ymysts = np.array([]); deltp_ps = np.array([]); sigmas = np.array([]); shots2 = np.array([])
plunges = np.array([]); rs = np.array([]); ls = np.array([])
shots = [190440, 190442, 190484, 190485, 190486, 184267, 167193, 167195, 184527, 187111]
for shot in shots:
    for plunge in [1, 2]:

        # Bad data.
        if shot == 187111 and plunge == 1:
            continue
        if shot == 167195 and plunge == 1:
            continue

        # Run normally, then pull out the solved vr and time-averaged blob vr. Grab some other things while we're
        # # here as well.
        print("#{} Plunge #{}".format(shot, plunge))
        bfs = BlobbyFarSOL.main(shot, plunge, showplot=False, min_npeaks=5)

        # Mask to restrict analysis to data within both RCP analysis domains.
        min_r = max(bfs.rcp_r.min(), bfs.blob_r.min())
        max_r = min(bfs.rcp_r.max(), bfs.blob_r.max())
        blob_mask = np.logical_and(bfs.blob_r > min_r, bfs.blob_r < max_r)

        # Extract some repeated variables for easier use.
        te = bfs.blob_te[blob_mask]
        ti = te * ti_mult
        ne = bfs.blob_ne[blob_mask]
        R = bfs.blob_r[blob_mask] / 100 # cm to m
        B = bfs.blob_b[blob_mask]
        vr = bfs.blob_vr[blob_mask]
        nn = bfs.neut_dens2[blob_mask]

        # Interpolation functions needed for later.
        f_lambda_ne = interp1d(bfs.rcp_r, bfs.lambda_ne)
        f_conn = interp1d(bfs.conns_r, bfs.conns_l_otf + bfs.conns_l_itf)
        f_conn_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
        f_rcp_ne = interp1d(bfs.rcp_r, bfs.rcp_ne)
        f_rcp_te = interp1d(bfs.rcp_r, bfs.rcp_te)

        # Calculate the normalized collisionality.
        vi = np.sqrt(3 * elec * te / md)  # m/s
        vte = np.sqrt(2 * elec * te / md)  # m/s
        nu_ei = (1 / (4 * np.pi) * ne * elec ** 4 / (eps ** 2 * me * md) * ln_lambda) * \
                (1 / (np.power(vi, 3) + 1.3 * np.power(vte, 3)))  # 1/s
        elec_gyrofreq = 1.609e-19 * B / me
        ion_gyrofreq = 1.609e-19 * B / md
        ion_gyrorad = np.sqrt(2 * 2.014 * 1.66e-27 * elec * te) / (elec * B)
        norm_coll = nu_ei * f_conn_otf(R * 100) / (elec_gyrofreq * ion_gyrorad)

        # Now calculate the normalized radius.
        ab = bfs.blob_drad[blob_mask] / 2 / 100  # cm to m
        Lotf = f_conn_otf(R * 100)
        L = f_conn(R * 100)
        # ahat = ab * np.power(R, 1 / 5) / (np.power(Lotf, 2 / 5) * np.power(ion_gyrorad, 4 / 5))
        ahat = ab * np.power(R, 1 / 5) / (np.power(L, 2 / 5) * np.power(ion_gyrorad, 4 / 5))
        norm_rad = np.power(ahat, 5 / 2)

        # Now the normalized velocity.
        cs1 = bfs.blob_cs[blob_mask]
        # norm_vr = vr / (np.power(2 * Lotf * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs1)
        norm_vr = vr / (np.power(2 * L * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs1)

        # Calculate the blob density and pressure above the background value.
        deltn_n = (ne - f_rcp_ne(R * 100)) / f_rcp_ne(R * 100)
        p = (te + ti) * ne
        rcp_p = (f_rcp_te(R * 100) + f_rcp_te(R * 100) * ti_mult) * f_rcp_ne(R * 100)
        deltp_p = (p - rcp_p) / rcp_p

        # Calculate the ionization and CX collision rates. Don't need the ionization one here but why not.
        nu_iz = nn * bfs.blob_ioniz_rate_coeffs[blob_mask]
        nu_cx = nn * bfs.blob_cx_rate_coeffs[blob_mask]

        # Now we should be ready to put it all together and calculate the predicted velocity from Tess's full
        # derivation.
        cs0 = np.sqrt(te / md_ev)
        rhos0 = cs0 / ion_gyrofreq
        cs2 = np.sqrt(1 + 3 * ti / te) * cs0

        def gen_vb(ab, sigma=2):
            yint = np.sqrt(2) * cs2 / np.sqrt(R * ab)
            instab_ratio = yint / yint
            numerator = np.sqrt(1 + ti / te) * np.sqrt(2 * ab / R) * cs0 * deltp_p
            denominator = instab_ratio + np.sqrt(2 * R * (1 + 3 * ti / te)) * sigma * np.power(ab, 5 / 2) / \
                          (np.sqrt(1 + ti / te) * L * np.square(rhos0)) + np.sqrt(R * ab / (2 * (1 + ti / te))) \
                          * nu_cx / cs0
            return numerator / denominator

        def calc_mystery_instab(sigma=2):
            numerator = np.sqrt(1 + ti / te) * np.sqrt(2 * ab / R) * cs0 * deltp_p
            instab_ratio = numerator / vr - np.sqrt(2 * R * (1 + 3 * ti / te)) * sigma * np.power(ab, 5 / 2) / \
                          (np.sqrt(1 + ti / te) * L * np.square(rhos0)) - np.sqrt(R * ab / (2 * (1 + ti / te))) \
                          * nu_cx / cs0
            yint = np.sqrt(2) * cs2 / np.sqrt(R * ab)
            return instab_ratio * yint

        def calc_sigma():
            numerator = np.sqrt(1 + ti / te) * np.sqrt(2 * ab / R) * cs0 * deltp_p
            tmp = numerator / vr - 1 - np.sqrt(R * ab / (2 * (1 + ti / te))) * nu_cx / cs0
            sigma = tmp / (np.sqrt(2 * R * (1 + 3 * ti / te)) * np.power(ab, 5 / 2) / \
                          (np.sqrt(1 + ti / te) * L * np.square(rhos0)))
            return sigma

        vr_pred = gen_vb(ab)
        # norm_vr_pred = vr_pred / (np.power(2 * Lotf * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs1)
        norm_vr_pred = vr_pred / (np.power(2 * L * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs1)

        # Put everything into our total arrays.
        vrs = np.append(vrs, vr)
        vrs_pred = np.append(vrs_pred, vr_pred)
        abs = np.append(abs, ab)
        norm_vrs = np.append(norm_vrs, norm_vr)
        norm_vrs_pred = np.append(norm_vrs_pred, norm_vr_pred)
        ahats = np.append(ahats, ahat)
        ymysts = np.append(ymysts, calc_mystery_instab())
        deltp_ps = np.append(deltp_ps, deltp_p)
        sigmas = np.append(sigmas, calc_sigma())
        shots2 = np.append(shots2, np.full(len(vr), shot))
        plunges = np.append(plunges, np.full(len(vr), plunge))
        rs = np.append(rs, R)
        ls = np.append(ls, L)

        # We can put aside a particular plunge so we can highlight on the final plot.
        # if shot == shot_plot and plunge == plunge_plot:
        #     plot_ab = ab
        #     plot_vr = vr
        #     plot_ab_pred = np.linspace(0, 0.03, 100)
        #     plot_vr_pred = gen_vb(plot_ab_pred)

        # Special case for the oft-studied 167195.
        if shot == 167195 and plunge == 2:
            gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
            with open(gfile_path, "rb") as f:
                gfile = pickle.load(f)
            R2 = gfile["R"]
            Z = gfile["Z"]
            Rs, Zs = np.meshgrid(R2, Z)
            psin = gfile["PSIRZ_NORM"]
            rcp_coord = zip(R, np.full(len(R), -0.185))
            rcp_psin_167195 = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))
            rcp_sigma_167195 = calc_sigma()

        elif shot == 190484:
            gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle"
            with open(gfile_path, "rb") as f:
                gfile = pickle.load(f)
            R2 = gfile["R"]
            Z = gfile["Z"]
            Rs, Zs = np.meshgrid(R2, Z)
            psin = gfile["PSIRZ_NORM"]
            rcp_coord = zip(R, np.full(len(R), -0.185))
            if plunge == 1:
                rcp_psin_190484_1 = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))
                rcp_sigma_190484_1 = calc_sigma()
            elif plunge == 2:
                rcp_psin_190484_2 = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))
                rcp_sigma_190484_2 = calc_sigma()

mask = ~np.isnan(ymysts)
p = np.polyfit(deltp_ps[mask], ymysts[mask]*1e-7, 1)
z = np.poly1d(p)
ymyst_fitx = np.linspace(0, 6, 2)
ymyst_fity = z(ymyst_fitx) * 1e7

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 8))

# First, a plot of the equation lines with my experimental lines.
ax1.scatter(abs, vrs, marker="^", s=75, edgecolor="k", zorder=15, color="tab:red")
# ax1.scatter(plot_ab, plot_vr, marker="^", s=75, edgecolor="k", zorder=15, color="tab:red")
# ax1.plot(plot_ab_pred, plot_vr_pred, color="tab:red")
ax1.grid(alpha=0.3, zorder=5)
ax1.set_xlabel(r"$\mathdefault{a_b}$ (m)", fontsize=14)
ax1.set_ylabel(r"$\mathdefault{v_r^{exp}}$ (m/s)", fontsize=14)

# Then a plot comparing the experimentally measured radial velocities to the predicted one from the equation.
ax2.plot([0, 2500], [0, 2500], color="k", linestyle="--")
sc = ax2.scatter(vrs_pred, vrs, c=ahats, marker="^", s=75, cmap="magma", edgecolor="k", zorder=15, vmin=0, vmax=3)
cbar = fig.colorbar(sc, ax=ax2)
ax2.grid(alpha=0.3, zorder=5)
ax2.set_xlabel(r"$\mathdefault{v_r^{pred}}$ (m/s)", fontsize=14)
ax2.set_ylabel(r"$\mathdefault{v_r^{exp}}$ (m/s)", fontsize=14)
cbar.set_label(r"$\mathdefault{\hat{a}}$", fontsize=14)

# Now a plot of the "mystery instability" vs. ahat since there seems to be a trend.
ax3.plot(ymyst_fitx, ymyst_fity, linestyle="--", color="k")
ax3.scatter(deltp_ps, ymysts, color="tab:red", edgecolors="k", marker="^", s=75, zorder=15)
ax3.grid(alpha=0.3, zorder=5)
ax3.set_xlabel(r"$\mathdefault{\delta_p \backslash \ p}$", fontsize=14)
ax3.set_ylabel(r"$\mathdefault{\gamma_{inst}}$ (s)", fontsize=14)
ax3.set_xlim([0, 6])
ax3.set_ylim([0, 1.5e7])

# Same as ax2, just colorbar is the needed sigma to force agreement.
ax4.plot([0, 2500], [0, 2500], color="k", linestyle="--")
sc = ax4.scatter(vrs_pred, vrs, c=sigmas, marker="^", s=75, cmap="magma", edgecolor="k", zorder=15, vmin=0, vmax=10)
cbar = fig.colorbar(sc, ax=ax4)
ax4.grid(alpha=0.3, zorder=5)
ax4.set_xlabel(r"$\mathdefault{v_r^{pred}}$ (m/s)", fontsize=14)
ax4.set_ylabel(r"$\mathdefault{v_r^{exp}}$ (m/s)", fontsize=14)
cbar.set_label(r"Needed $\mathdefault{\sigma}$ to force agreement", fontsize=14)

fig.tight_layout()
fig.show()
