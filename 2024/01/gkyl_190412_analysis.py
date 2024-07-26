import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg
import pandas as pd
from tqdm import tqdm

# Any constants up front.
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
        # mp = 1.672623e-27
        # AMU = 2.014  # Deuterium ions
        # mi = mp * AMU
        # me = 9.10938188e-31
        eV = 1.602e-19

        # c_s = np.sqrt(40 * eV / mi)

        # For equilibrium profiles
        d = pg.GData("%s_electron_M0_%d.bp" % (fp, fstart))
        dg = pg.GInterpModal(d, 1, 'ms')
        X, nElc = dg.interpolate(0)
        # nElc = nElc[:, :, :, 0]

        # zi = len(X[2]) // 2  # to take a slice at midplane

        # Nx = len(X[0]) - 1
        # Nz = len(X[2]) - 1
        # Ny = len(X[1]) - 1

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
            if t == fend - 1:
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

            # Phi eq profile calculation
            if not skip_phi:
                d = pg.GData("%s_phi_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, phi = dg.interpolate(0)
                # phiZmin = phi[:, :, 0, 0]
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
        def bmag(x_val):
            B_axis = 2.04
            R = 2.30
            R0 = 1.722
            B0 = B_axis * (R0 / R)
            return B0 * R / x_val

        BB = np.zeros(self.ne.shape[1:])
        for i in range(0, len(self.r)):
            BB[i, :, :] = bmag(self.r[i])
        self.BB = BB

        print("Calculating the magnetic field gradient...")
        gradB = np.gradient(BB, x, y, z)
        self.gradBx = gradB[0]  # Units of T / m
        self.gradBy = gradB[1]  # ''
        self.gradBz = gradB[2]  # Units of T / radian


# Load plasma into object.
plasma = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-190412-v2/", "d3d-190412-v2", 600, 999,
                 skip_ti=False, skip_phi=False, skip_upar=True)
midp = len(plasma.p) // 2
midz = len(plasma.z) // 2

# Load RCP data.
rcp_rsep = 2.23119  # From EFITviewer
rcp_path1 = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP190412_1.tab"
rcp1 = pd.read_csv(rcp_path1, delimiter="\t")
rcp_r1 = rcp1["R(cm)"].values / 100
rcp_te1 = rcp1["Te(eV)"].values
rcp_ne1 = rcp1["Ne(E18 m-3)"].values * 1e18
rcp_path2 = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP190412_2.tab"
rcp2 = pd.read_csv(rcp_path2, delimiter="\t")
rcp_r2 = rcp2["R(cm)"].values / 100
rcp_te2 = rcp2["Te(eV)"].values
rcp_ne2 = rcp2["Ne(E18 m-3)"].values * 1e18

# The R position of the separatrix at the Gkeyll midz. This is essentially Z=0 since that is where halfway between
# the upper and lower baffle occurs. Room to play with this since the simple helical approximation makes this tough.
gk_rsep = 2.241

# Average density over the simulation time.
gk_r = plasma.r
gk_z = plasma.z
gk_ne_avg = plasma.ne[:, :, midp, midz].mean(axis=0)
gk_te_avg = plasma.te[:, :, midp, midz].mean(axis=0)
gk_dndr_avg = np.gradient(gk_ne_avg, gk_r)

end_idx1 = 20
end_idx2 = 10
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6))
# ax1.plot(rcp_r1[end_idx1:]-rcp_rsep, rcp_ne1[end_idx1:], label="RCP", color="tab:red", lw=2)
ax1.plot(rcp_r2[end_idx2:]-rcp_rsep, rcp_ne2[end_idx2:], color="tab:red", lw=2)
ax1.plot(gk_r[8:-2]-gk_rsep, gk_ne_avg[8:-2], color="k", lw=3)
ax1.plot(gk_r[8:-2]-gk_rsep, gk_ne_avg[8:-2], label="Gkeyll", color="tab:purple", lw=2)
ax1.legend()
# ax1.set_xlabel("R (m)", fontsize=12)
ax1.set_xticklabels([])
# ax1.set_xlim([0.03, 0.12])
ax1.set_ylabel(r"$\mathdefault{n_e}$ ($\mathdefault{m^{-3}}$)", fontsize=12)
# ax1.set_ylim([1e18, 5e19])
ax1.set_yscale("log")
ax1.grid(alpha=0.3, which="both")
# ax1.set_xlim(np.array([None, 2.40])
# ax1.set_xticks([2.25, 2.30, 2.35, 2.40])
# ax2.plot(rcp_r1[end_idx1:]-rcp_rsep, rcp_te1[end_idx1:], label="RCP", color="tab:red", lw=2)
ax2.plot(rcp_r2[end_idx2:]-rcp_rsep, rcp_te2[end_idx2:], color="tab:red", lw=2)
ax2.plot(gk_r-gk_rsep, gk_te_avg, color="k", lw=3)
ax2.plot(gk_r-gk_rsep, gk_te_avg, label="Gkeyll", color="tab:purple", lw=2)
ax2.set_xlabel(r"R-$\mathdefault{R_{sep}}$ (m)", fontsize=12)
ax2.set_ylabel(r"$\mathdefault{T_e}$ (eV)", fontsize=12)
ax2.set_yscale("log")
# ax1.set_xticklabels(["{:.2f}".format(n) for n in [2.25, 2.30, 2.35, 2.40]])
# ax2.set_xlim([0.03, 0.12])
ax2.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.show()
