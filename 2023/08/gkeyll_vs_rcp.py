# Script to see how our Gkeyll run for 167196 compares to the RCP.
import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg
import pandas as pd
from tqdm import tqdm

mD = 931.49e6 / 3e8 ** 2

class Plasma3:

    def __init__(self, root, fname, fstart, fend, skip_ti=False, skip_phi=False, skip_upar=False):
        """
        root (str): Path to the directory containing the background data.
        fstart (int): Starting frame number.
        fend (end): Ending frame number.
        """

        # Load in the plasma background data. This is taken from Tess's nstx-diagnostic.py script.
        step = 1
        Nt = fend - fstart
        Nt //= step
        dt = 5e-7

        fp = "{:}{:}".format(root, fname)
        print("Data location: {:}".format(fp))

        # physical constants
        mp = 1.672623e-27
        AMU = 2.014  # Deuterium ions
        mi = mp * AMU
        me = 9.10938188e-31
        eV = 1.602e-19

        c_s = np.sqrt(40 * eV / mi)

        # For equilibrium profiles
        d = pg.GData("%s_electron_M0_%d.bp" % (fp, fstart))
        dg = pg.GInterpModal(d, 1, 'ms')
        X, nElc = dg.interpolate(0)
        nElc = nElc[:, :, :, 0]

        zi = len(X[2]) // 2  # to take a slice at midplane

        Nx = len(X[0]) - 1
        Nz = len(X[2]) - 1
        Ny = len(X[1]) - 1

        # calculate grids
        x = X[0]
        dx = np.diff(x)[0]
        x = x + dx / 2
        x = x[:-1]

        y = X[1]
        dy = np.diff(y)[0]
        y = y + dy / 2
        y = y[:-1]

        z = X[2]
        dz = np.diff(z)[0]
        z = z + dz / 2
        z = z[:-1]

        # Create arrays for equilibrium profiles that are averaged in y, z and time
        nElc_tot = []
        tElc_tot = []
        tIon_tot = []
        upar_tot = []
        phi_tot = []
        elecx_tot = []
        elecy_tot = []
        elecz_tot = []
        print("Loading background plasma...")
        for t in range(fstart, fend, step):
            if t % 100 == 0:
                print('frame = ', t)
            if (t == fend - 1):
                print("frame = ", t)

            # electron density
            d = pg.GData("%s_electron_M0_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, nElc = dg.interpolate(0)
            nElc = nElc[:, :, :, 0]
            nElc_tot.append(nElc)

            # electron temperature
            d = pg.GData("%s_electron_Temp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tElc = dg.interpolate(0)
            tElc = tElc[:, :, :, 0] / eV
            tElc_tot.append(tElc)

            # ion temperature
            if not skip_ti:
                d = pg.GData("%s_ion_Temp_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, tIon = dg.interpolate(0)
                tIon = tIon[:, :, :, 0] / eV
                tIon_tot.append(tIon)

            # parallel flow velocity
            if not skip_upar:
                d = pg.GData("%s_ion_Upar_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, upar = dg.interpolate(0)
                upar = upar[:, :, :, 0]
                upar_tot.append(upar)

            # Phi eq profile calcuation
            if not skip_phi:
                d = pg.GData("%s_phi_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, phi = dg.interpolate(0)
                phiZmin = phi[:, :, 0, 0]
                phi = phi[:, :, :, 0]
                phi_tot.append(phi)

                # The electric field as the gradient of the potential.
                elec = np.gradient(-phi, x, y, z)
                elecx_tot.append(elec[0])
                elecy_tot.append(elec[1])
                elecz_tot.append(elec[2])

        nElc_tot = np.array(nElc_tot)
        tElc_tot = np.array(tElc_tot)
        tIon_tot = np.array(tIon_tot)
        upar_tot = np.array(upar_tot)
        phi_tot = np.array(phi_tot)
        elecx_tot = np.array(elecx_tot)
        elecy_tot = np.array(elecy_tot)
        elecz_tot = np.array(elecz_tot)
        self.ne = nElc_tot
        self.te = tElc_tot
        self.ti = tIon_tot
        self.upar = upar_tot
        self.vp = phi_tot
        self.er = elecx_tot
        self.ep = elecy_tot
        self.ez = elecz_tot
        self.r = x
        self.p = y
        self.z = z
        self.dt = dt
        if not skip_ti:
            self.cs = np.sqrt((self.te + self.ti) / mD)

        # Calculate the B field and its gradients.
        def bmag(x):
            B_axis = 2.04
            R = 2.30
            R0 = 1.722
            B0 = B_axis * (R0 / R)
            return B0 * R / x
        BB = np.zeros(self.ne.shape[1:])
        for i in range(0, len(self.r)):
            BB[i, :, :] = bmag(self.r[i])
        self.BB = BB

        print("Calculating the magnetic field gradient...")
        gradB = np.gradient(BB, x, y, z)
        self.gradBx = gradB[0]  # Units of T / m
        self.gradBy = gradB[1]  # ''
        self.gradBz = gradB[2]  # Units of T / radian


plasma = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-167196-v4/", "d3d-167196-v4", 600, 999,  # Temporary, go back to 999
                 skip_ti=False, skip_phi=False)  # 0-999
midp = len(plasma.p) // 2
midz = len(plasma.z) // 2

# The R position of the separatrix at the Gkeyll midz (Z = 0 essentially).
gk_rsep = 2.259

# Average density over the simulation time.
gk_r = plasma.r
gk_z = plasma.z
gk_ne_avg = plasma.ne[:, :, midp, midz].mean(axis=0)
gk_te_avg = plasma.te[:, :, midp, midz].mean(axis=0)
gk_dndr_avg = np.gradient(gk_ne_avg, gk_r)

# The parallel flux.
gk_upar_avg = plasma.upar[:, :, midp, midz].mean(axis=0)
gk_gamma_par = gk_ne_avg * gk_upar_avg
gk_dgdr_par = np.gradient(gk_gamma_par, gk_r)
gk_dgdz_par = np.gradient(plasma.ne[:, :, midp, :].mean(axis=0) * plasma.upar[:, :, midp, :].mean(axis=0), gk_z, axis=1)[:, midz]
gk_dgdr_rad = -gk_dgdr_par - gk_dgdz_par
gk_gamma_rad = np.trapz(gk_dgdr_rad, gk_r)

# Estimate for radial diffusion coefficient (assuming Gamma_par = Gamma_r)
gk_deff = -gk_dndr_avg / gk_gamma_rad
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.plot(gk_r, gk_deff)
ax1.set_yscale("log")
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Deff (m2/s)")
ax1.grid(alpha=0.3)
fig.tight_layout()
fig.show()

# Plot considering a synthetic probe measuring blobs.
sig_locr = np.linspace(2.32, 2.38, 5)
sig_locp = np.zeros(len(sig_locr))
sig_t = np.arange(0, plasma.ne.shape[0] * plasma.dt, plasma.dt)
sig_ne = np.zeros((sig_locr.shape[0], plasma.ne.shape[0]))
sig_ep = np.zeros((sig_locr.shape[0], plasma.ne.shape[0]))
for f in range(0, plasma.ne.shape[0]):
    ne = plasma.ne[f]
    # ep = plasma.ep[f]
    for i in range(0, len(sig_locr)):
        r = sig_locr[i]
        # p = sig_locp[i]
        idx_r = np.argmin(np.abs(r - plasma.r))
        # idx_p = np.argmin(np.abs(p - plasma.p))
        sig_ne[i, f] = ne[idx_r, midp, midz]
        # sig_ep[i, f] = ne[idx_r, midp, midz]
rms_ne = np.array([np.sqrt(np.mean(np.square(sig_ne[i]))) for i in range(0, sig_ne.shape[0])])

# A super-duper rough estimate of blob counting.
blob_counts = np.zeros(rms_ne.shape)
tmp_blob_ne = np.zeros(rms_ne.shape)
tmp_counts = np.zeros(rms_ne.shape)
for i in range(0, sig_ne.shape[0]):
    ne = sig_ne[i]
    in_peak = False
    for j in range(1, len(ne)):

        # Criteria for start of a peak.
        if ne[j-1] < rms_ne[i] and ne[j] > rms_ne[i]:
            in_peak = True
            blob_counts[i] += 1
            tmp_blob_ne[i] += ne[j]
            tmp_counts[i] += 1

        # Criteria for in a peak and still in the peak. Do not count blob_counts again, would be doable-counting the
        # same peak. blob_ne gets added since we use tmp_counts to get an average of all the data point included at the
        # end.
        if ne[j-1] > rms_ne[i] and ne[j] > rms_ne[i]:
            in_peak = True
            tmp_blob_ne[i] += ne[j]
            tmp_counts[i] += 1

        # Criteria for end of a peak.
        if ne[j-1] > rms_ne[i] and ne[j] < rms_ne[i]:
            in_peak = False
blob_ne = tmp_blob_ne / tmp_counts
blob_fb = blob_counts / sig_t[-1]

# RCP data from 167195. Data beyond 233 should be removed.
rcp_shift = 0.0
print("IMPORTANT: Dmitry has said suggested eeither a 1 or 2 cm inward shift" +\
        " of the RCP data to bring it into better alignment with TS. The " + \
        "current shift is:")
print("Shift = {:.2f} cm".format(rcp_shift))
rcp_rsep = 2.216
rcp1 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP167195_1.tab", delimiter="\t")
rcp2 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP167195_2.tab", delimiter="\t")
rcp1 = rcp1[rcp1["R(cm)"] < 233]
rcp2 = rcp2[rcp2["R(cm)"] < 233]
rcp_r1 = rcp1["R(cm)"] / 100 + rcp_shift
rcp_r2 = rcp2["R(cm)"] / 100 + rcp_shift
rcp_ne1 = rcp1["Ne(E18 m-3)"] * 1e18
rcp_ne2 = rcp2["Ne(E18 m-3)"] * 1e18
rcp_te1 = rcp1["Te(eV)"]
rcp_te2 = rcp2["Te(eV)"]

# Conditional averaging data as well.
ca1 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_167195_1.tab", delimiter="\t")
ca2 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_167195_2.tab", delimiter="\t")
ca_r1 = ca1["R(cm)"] / 100 + rcp_shift
ca_r2 = ca2["R(cm)"] / 100 + rcp_shift
ca_dt1 = (ca1["Time(ms)"].iloc[1] - ca1["Time(ms)"].iloc[0]) / 1000  # ms to s
ca_dt2 = (ca2["Time(ms)"].iloc[1] - ca2["Time(ms)"].iloc[0]) / 1000
ca_fb1 = ca1["Npeaks"] / ca_dt1
ca_fb2 = ca2["Npeaks"] / ca_dt2

colors = ["C{}".format(i) for i in range(0, 10)]
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 7))
ax1.plot(rcp_r1-rcp_rsep, rcp_ne1, label="RCP #1")
ax1.plot(rcp_r2-rcp_rsep, rcp_ne2, label="RCP #2")
ax1.plot(gk_r-gk_rsep, gk_ne_avg, label="Gkeyll")
ax1.legend()
ax1.set_xlabel("R (m)")
ax1.set_ylabel(r"$\mathdefault{n_e}$ (m-3)")
ax1.set_ylim([5e17, 1e20])
ax1.set_yscale("log")
ax1.grid(alpha=0.3, which="both")
ax2.plot(rcp_r1-rcp_rsep, rcp_te1, label="RCP #1")
ax2.plot(rcp_r2-rcp_rsep, rcp_te2, label="RCP #2")
ax2.plot(gk_r-gk_rsep, gk_te_avg, label="Gkeyll")
ax2.set_xlabel(r"R-$\mathdefault{R_{sep}}$ (m)")
ax2.set_ylabel(r"$\mathdefault{T_e}$ (eV)")
ax2.set_yscale("log")
ax2.grid(alpha=0.3, which="both")
for i in range(0, sig_ne.shape[0]):
    ax3.plot(sig_t * 1e6, sig_ne[i], label="{:.3f}".format(sig_locr[i]), color=colors[i])
    ax3.axhline(rms_ne[i], linestyle="--", color=colors[i])
ax3.set_xlabel("Time (us)")
ax3.legend()
ax4.plot(ca_r1, ca_fb1)
ax4.plot(ca_r2, ca_fb2)
ax4.plot(sig_locr, blob_fb)
ax4.set_yscale("log")
ax4.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.show()

end_idx1 = -6
end_idx2 = -6
# fontsize = 12  # For paper
fontsize = 16  # For presentations
#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
ax1.plot(rcp_r1[:end_idx1]-rcp_rsep, rcp_ne1[:end_idx1], label="RCP", color="tab:red", lw=2)
ax1.plot(rcp_r2[:end_idx2]-rcp_rsep, rcp_ne2[:end_idx2], color="tab:red", lw=2)
ax1.plot(gk_r[8:-2]-gk_rsep, gk_ne_avg[8:-2], color="k", lw=3)
ax1.plot(gk_r[8:-2]-gk_rsep, gk_ne_avg[8:-2], label="Gkeyll", color="tab:purple", lw=2)
ax1.legend()
ax1.set_xlabel(r"R-$\mathdefault{R_{sep}}$ (m)", fontsize=fontsize)
# ax1.set_xticklabels([])
ax1.set_xlim([0.03, 0.12])
ax1.set_ylabel(r"$\mathdefault{n_e}$ ($\mathdefault{m^{-3}}$)", fontsize=fontsize)
ax1.set_ylim([1e18, 5e19])
ax1.set_yscale("log")
ax1.grid(alpha=0.3, which="both")
# ax1.set_xlim(np.array([None, 2.40])
# ax1.set_xticks([2.25, 2.30, 2.35, 2.40])
ax2.plot(rcp_r1[:end_idx1]-rcp_rsep, rcp_te1[:end_idx1], label="RCP", color="tab:red", lw=2)
ax2.plot(rcp_r2[:end_idx2]-rcp_rsep, rcp_te2[:end_idx2], color="tab:red", lw=2)
ax2.plot(gk_r-gk_rsep, gk_te_avg, color="k", lw=3)
ax2.plot(gk_r-gk_rsep, gk_te_avg, label="Gkeyll", color="tab:purple", lw=2)
ax2.set_xlabel(r"R-$\mathdefault{R_{sep}}$ (m)", fontsize=fontsize)
ax2.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=fontsize)
ax2.set_yscale("log")
# ax1.set_xticklabels(["{:.2f}".format(n) for n in [2.25, 2.30, 2.35, 2.40]])
ax2.set_xlim([0.03, 0.12])
ax2.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.show()

# Histogram of the fluctuations at each x location.
bins = np.linspace(-7000, 7000, 15)
all_counts = {}
avgs = []
stds = []
for ir in range(0, plasma.ep.shape[1]):
    counts = np.zeros(len(bins) - 1)
    for t in range(0, plasma.ep.shape[0]):
        ep = plasma.ep[t, :, :, midz]
        tmp_counts, _ = np.histogram(ep[ir, :], bins=bins)
        counts += tmp_counts
    all_counts[ir] = counts
    avgs.append(plasma.ep[:, ir, :, midz].mean())
    stds.append(plasma.ep[:, ir, :, midz].std())
avgs = np.array(avgs)
stds = np.array(stds)

barx = [(bins[i+1] + bins[i]) / 2 for i in range(0, len(bins)-1)]
fig, ax1 = plt.subplots()
ax1.fill_between(plasma.r, avgs-stds, avgs+stds, alpha=0.3)
ax1.plot(plasma.r, avgs)
fig.tight_layout()
fig.show()

# Histogram of the radial velocity values at each x location.
bins_dict = {}
imp_Z = 5
imp_mass = 184.84 * 1.66E-27  # amu x kg
e = 1.609e-19
conn_length = 10
lamb_rs = np.zeros(len(plasma.r) - 1)
tidx = 10
for ir in range(0, len(plasma.r)):

    # Extract needed variables...
    ep = plasma.ep[:, ir, :, midz]
    er = plasma.er[:, ir, :, midz]
    b = plasma.BB[ir, midp, midz]
    omega_c = imp_Z * e * b / imp_mass
    cs = plasma.cs[:, ir, :, midz]

    # Time derivative of Er for the polarization drift.
    der_dt = np.zeros(er.shape)
    for i in range(0, len(der_dt)):
        if i == 0:
            der_dt[i] = (er[i+1] - er[i]) / plasma.dt
        elif i == len(der_dt) - 1:
            der_dt[i] = (er[i] - er[i-1]) / plasma.dt
        else:
            der_dt[i] = (er[i+1] - er[i-1]) / (2 * plasma.dt)

    # Calculate the radial velocity.
    vr = ep / b + 1 / (omega_c * b) * der_dt

    # Bin the data.
    bins = np.linspace(vr.min()*0.75, vr.max()*0.75, 21)
    hist, bin_edges = np.histogram(vr, bins)
    bins_dict[ir] = {"hist":hist, "bin_edges":bin_edges}

    # Now go through one radial location at a time, calculating lambda for W5+ ions. First calculate the critical vr
    # for an ion to make it to the next cell.
    if ir == len(plasma.r) - 1:
        continue
    delt_r = plasma.r[ir+1] - plasma.r[ir]
    vr_crit = delt_r * cs.mean() / conn_length

    # Then see what the probability of a velocity greater than or equal to it is. Use the closest bin as an
    # approximation.
    # idx = np.argmin(np.abs(bin_edges - vr_crit))
    # p_star = hist[idx:].sum() / hist.sum()

    p_star = np.count_nonzero(vr > vr_crit) / np.count_nonzero(vr)

    # Then calculate lambda.
    lamb_r = delt_r / (1 - p_star)
    lamb_rs[ir] = lamb_r

fig, ax1 = plt.subplots()
ax1.plot(plasma.r[:-1], lamb_rs)
ax1.set_xlabel("R (m)")
ax1.set_ylabel("lambda_r (m)")
fig.tight_layout()
fig.show()
