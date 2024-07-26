# Script to see how our Gkeyll run for 167196 compares to the RCP.
import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg
import pandas as pd
from tqdm import tqdm
from scipy.stats import skew, kurtosis

mD = 931.49e6 / 3e8 ** 2


class Plasma3:

    def __init__(self, root, fname, fstart, fend, skip_ti=False, skip_phi=False, skip_upar=False, tung_flag=False):
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
        nCar_tot = []
        tCar_tot = []
        nIon_tot = []
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

            # If a carbon run load additional carbon files.
            if tung_flag:

                # carbon density
                d = pg.GData("%s_tungsten_M0_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, nCar = dg.interpolate(0)
                nCar = nCar[:, :, :, 0]
                nCar_tot.append(nCar)

                # carbon temperature
                d = pg.GData("%s_tungsten_Temp_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, tCar = dg.interpolate(0)
                tCar = tCar[:, :, :, 0] / eV
                tCar_tot.append(tCar)

                # and the main ion density since ne=ni is no longer true.
                d = pg.GData("%s_ion_M0_%d.bp" % (fp, t))
                dg = pg.GInterpModal(d, 1, 'ms')
                X, nIon = dg.interpolate(0)
                nIon = nIon[:, :, :, 0]
                nIon_tot.append(nIon)

        nElc_tot = np.array(nElc_tot)
        tElc_tot = np.array(tElc_tot)
        nIon_tot = np.array(nIon_tot)
        tIon_tot = np.array(tIon_tot)
        nCar_tot = np.array(nCar_tot)
        tCar_tot = np.array(tCar_tot)
        upar_tot = np.array(upar_tot)
        phi_tot = np.array(phi_tot)
        elecx_tot = np.array(elecx_tot)
        elecy_tot = np.array(elecy_tot)
        elecz_tot = np.array(elecz_tot)
        self.ne = nElc_tot
        self.ni = nIon_tot
        self.te = tElc_tot
        self.ti = tIon_tot
        self.nc = nCar_tot  # Lazy, should rename all the Car's above
        self.tc = tCar_tot
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


# Gkyell run with only e and D, matched to RCP data pretty well.
plasma = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-167196-v4/", "d3d-167196-v4", 600, 999,
                 skip_ti=False, skip_phi=False, skip_upar=True)  # 0-999

# Gkyell run where I have added a 2% C3+ background to the above run.
plasma_c3 = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-167196-v4-w5/", "d3d-167196-v4-w5", 600, 994,
                    skip_ti=False, skip_phi=False, skip_upar=True, tung_flag=True)  # 0-999

# Same grid for both, same midpoints.
midp = len(plasma.p) // 2
midz = len(plasma.z) // 2

# Density of each species and it's relative fluctuation levels.
ne_avg = plasma.ne[:, :, midp, midz].mean(axis=0)  # First mean poloidal, second is time
ne_rms = np.sqrt(np.square(plasma.ne[:, :, midp, midz] - ne_avg).mean(axis=0))
ne_avg_c3 = plasma_c3.ne[:, :, midp, midz].mean(axis=0)  # First mean poloidal, second is time
ne_rms_c3 = np.sqrt(np.square(plasma_c3.ne[:, :, midp, midz] - ne_avg_c3).mean(axis=0))
ni_avg_c3 = plasma_c3.ni[:, :, midp, midz].mean(axis=0)
nc_avg_c3 = plasma_c3.nc[:, :, midp, midz].mean(axis=0)

# Temperatures
te_avg = plasma.te[:, :, midp, midz].mean(axis=0)
ti_avg = plasma.ti[:, :, midp, midz].mean(axis=0)
te_avg_c3 = plasma_c3.te[:, :, midp, midz].mean(axis=0)
ti_avg_c3 = plasma_c3.ti[:, :, midp, midz].mean(axis=0)
tc_avg_c3 = plasma_c3.tc[:, :, midp, midz].mean(axis=0)

# Skewness and excess kurtosis.
skewness = skew(np.sqrt(np.square(plasma.ne[:, :, midp, midz] - ne_avg) - ne_avg))
skewness_c3 = skew(np.sqrt(np.square(plasma_c3.ne[:, :, midp, midz] - ne_avg_c3) - ne_avg_c3))
kurt = kurtosis(np.sqrt(np.square(plasma.ne[:, :, midp, midz] - ne_avg) - ne_avg))
kurt_c3 = kurtosis(np.sqrt(np.square(plasma_c3.ne[:, :, midp, midz] - ne_avg_c3) - ne_avg_c3))

std_color = "tab:purple"
car_color = "tab:red"
fig, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2, 3, sharex=True, figsize=(10, 5))

ax1.plot(plasma.r, (ne_avg - ne_rms) / ne_avg, label="Standard", color=std_color)
ax1.plot(plasma.r, (ne_avg_c3 - ne_rms_c3) / ne_avg_c3, label="2% C3+", color=car_color)
ax1.set_ylabel("(n - <n>) / <n>")
ax1.legend()

ax2.plot(plasma.r, te_avg, color=std_color)
ax2.plot(plasma.r, ti_avg, color=std_color, linestyle="--")
ax2.plot(plasma.r, te_avg_c3, color=car_color, label="Te")
ax2.plot(plasma.r, ti_avg_c3, color=car_color, label="Ti", linestyle="--")
ax2.plot(plasma.r, tc_avg_c3, color=car_color, label="Tc", linestyle=":")
ax2.legend()
ax2.set_ylabel("Temperature (eV)")

ax3.plot(plasma.r, ne_avg, color=std_color)
ax3.plot(plasma.r, ne_avg_c3, color=car_color, label="ne")
ax3.plot(plasma.r, ni_avg_c3, color=car_color, linestyle="--", label="ni")
ax3.plot(plasma.r, nc_avg_c3, color=car_color, linestyle=":", label="nc")
ax3.set_ylabel("Density (m-3)")
ax3.legend()

ax4.plot(plasma.r, nc_avg_c3 / ne_avg_c3 * 100, color=car_color)
ax4.set_ylabel("nC / nD (%)")

ax5.plot(plasma.r, skewness, color=std_color)
ax5.plot(plasma.r, skewness_c3, color=car_color)
ax5.set_ylabel("Skewness")

ax6.plot(plasma.r, kurt, color=std_color)
ax6.plot(plasma.r, kurt_c3, color=car_color)
ax6.set_ylabel("Kurtosis")

fig.tight_layout()
fig.show()
