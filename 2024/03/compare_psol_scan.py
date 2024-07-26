import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg
import pandas as pd
from tqdm import tqdm

mD = 931.49e6 / 3e8 ** 2

class Plasma3:

    def __init__(self, root, fname, fstart, fend, skip_ti=False, skip_phi=False, skip_upar=False, skip_te=False, skip_ne=False):
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
            if not skip_ne:
                d = pg.GData("%s_electron_M0_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, nElc = dg.interpolate(0)
                nElc = nElc[:, :, :, 0]
                nElc_tot.append(nElc)

            # electron temperature
            if not skip_te:
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
        if not skip_ti or not skip_te:
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


# Load each of the runs. 
pscan1 = Plasma3("/Users/zamperini/gkyldir/d3d-167196-v4-psol-scan1/", "d3d-167196-v4-psol-scan1", 600, 999,  
                 skip_ti=True, skip_phi=True, skip_te=True, skip_upar=True)  
pscan2 = Plasma3("/Users/zamperini/gkyldir/d3d-167196-v4-psol-scan2/", "d3d-167196-v4-psol-scan2", 600, 999,  
                 skip_ti=True, skip_phi=True, skip_te=True, skip_upar=True)
pscan3 = Plasma3("/Users/zamperini/gkyldir/d3d-167196-v4-psol-scan3/", "d3d-167196-v4-psol-scan3", 600, 999,
                 skip_ti=True, skip_phi=True, skip_te=True, skip_upar=True)
pscan4 = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-167196-v4/", "d3d-167196-v4", 600, 999, 
                 skip_ti=True, skip_phi=True, skip_te=True, skip_upar=True)
midp = len(pscan1.p) // 2
midz = len(pscan1.z) // 2

def get_avg_ne(plasma):
    gk_r = plasma.r
    gk_ne_avg = plasma.ne[:, :, midp, midz].mean(axis=0)
    return gk_r, gk_ne_avg
    
r, ne1 = get_avg_ne(pscan1) 
r, ne2 = get_avg_ne(pscan2) 
r, ne3 = get_avg_ne(pscan3) 
r, ne4 = get_avg_ne(pscan4)    
    
fig, ax1 = plt.subplots(figsize=(5,4))
#ax1.plot(r, ne1, label="0.05")
ax1.plot(r, ne2, label="0.25")
ax1.plot(r, ne3, label="0.75")
ax1.plot(r, ne4, label="0.45")
ax1.set_yscale("log")
ax1.legend()
ax1.set_ylim([8e17, 5e19])
ax1.grid(alpha=0.25, which="both")
ax1.set_xlabel("R (m)", fontsize=14)
ax1.set_ylabel("ne (m-3)", fontsize=14)
fig.tight_layout()
fig.show()
    
    
