import numpy as np
import postgkyl as pg
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker, cm
from tqdm import tqdm
import scipy.integrate as integrate
import sys
from matplotlib.collections import LineCollection
from scipy.signal import savgol_filter
import pandas as pd
import time

mD = 931.49e6 / 3e8 ** 2


class Impurity:

    def __init__(self, charge, mass, r, p, z, dt, fbirth):
        """
        charge (int): charge of ion
        mass (float): mass of ion in kg
        r (float): R location (m)
        p (float): Binormal (~poloidal) location (m)
        """

        self.charge = charge
        self.mass = mass
        self.r = r
        self.p = p
        self.z = z
        self.vr = 0.0
        self.vp = 0.0
        self.vz = 0.0
        self.t = 0.0
        self.dt = dt
        self.rhist = [r]
        self.phist = [p]
        self.zhist = [z]
        self.vrhist = [0.0]
        self.vphist = [0.0]
        self.vzhist = [0.0]
        self.ephist = [0.0]
        self.killed = False
        self.fbirth = fbirth

    def calc_force(self):
        pass


class Plasma:

    def __init__(self, root, fstart, fend):
        """
        root (str): Path to the directory containing the background data.
        fstart (int): Starting frame number.
        fend (end): Ending frame number.
        """

        # Load in the plasma background data. This is taken from Tess's nstx-diagnostic.py script.
        step = 1
        Nt = fend - fstart
        Nt //= step
        dt = 1e-6

        fp = "{:}nstx-neut".format(root)
        print("Data location: {:}".format(fp))

        # physical constants
        mp = 1.672623e-27
        AMU = 2.014  # Deuterium ions
        mi = mp * AMU
        me = 9.10938188e-31
        eV = 1.602e-19

        c_s = np.sqrt(40 * eV / mi)

        # For equilibrium profiles
        d = pg.GData("%s_electron_GkM0_%d.bp" % (fp, fstart))
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
            d = pg.GData("%s_electron_GkM0_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, nElc = dg.interpolate(0)
            nElc = nElc[:, :, :, 0]
            nElc_tot.append(nElc)

            # electron temperature
            d = pg.GData("%s_electron_GkTemp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tElc = dg.interpolate(0)
            tElc = tElc[:, :, :, 0] / eV
            tElc_tot.append(tElc)

            # ion temperature
            d = pg.GData("%s_ion_GkTemp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tIon = dg.interpolate(0)
            tIon = tIon[:, :, :, 0] / eV
            tIon_tot.append(tIon)

            # Phi eq profile calcuation
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
        phi_tot = np.array(phi_tot)
        elecx_tot = np.array(elecx_tot)
        elecy_tot = np.array(elecy_tot)
        elecz_tot = np.array(elecz_tot)
        self.ne = nElc_tot
        self.te = tElc_tot
        self.ti = tIon_tot
        self.vp = phi_tot
        self.er = elecx_tot
        self.ep = elecy_tot
        self.r = x
        self.p = y
        self.z = z
        self.dt = dt
        self.cs = np.sqrt((self.te + self.ti) / mD)


class Plasma2:

    def __init__(self, root="/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-posD04-k135-3x2v/", tstart=0, tend=654):
        """
        Loads the plasma background from the d3d-posD04-k135-3xv2 case.
        """

        # These files have shape (Nt, Nx, Ny, Nz) = (653, 192, 192, 32)
        print("Loading...")
        print("  - ne")
        nElc = np.load(root + 'd3d_elc_density.npy')
        print("  - Te")
        tElc = np.load(root + 'd3d_elc_temp.npy')
        print("  - Ti")
        tIon = np.load(root + 'd3d_ion_temp.npy')
        print("  - phi")
        phi = np.load(root + 'd3d_phi.npy')

        # Time array
        tGrid = np.linspace(0.0, 652e-6, 653)

        # Grid in field-line following coordinate system.
        # Separatrix is at X = 0.1 m
        xGrid = np.linspace(0., 0.15, 192)
        yGrid = np.linspace(-7.074822e-02, 7.074822e-02, 192)
        zGrid = np.linspace(-np.pi, np.pi, 32)

        # The electric field as the gradient of the potential.
        print("Calculating electric field...")
        elecx_tot = np.zeros(phi.shape)
        elecy_tot = np.zeros(phi.shape)
        elecz_tot = np.zeros(phi.shape)
        for t in range(0, phi.shape[0]):
            elec = np.gradient(-phi[t], xGrid, yGrid, zGrid)
            elecx_tot[t] = elec[0]
            elecy_tot[t] = elec[1]
            elecz_tot[t] = elec[2]

        self.dt = tGrid[1] - tGrid[0]
        self.ne = nElc[tstart:tend]
        self.te = tElc[tstart:tend]
        self.ti = tIon[tstart:tend]
        self.vp = phi[tstart:tend]
        self.er = elecx_tot[tstart:tend]
        self.ep = elecy_tot[tstart:tend]
        self.ez = elecz_tot[tstart:tend]
        self.r = xGrid
        self.p = yGrid
        self.z = zGrid
        self.cs = np.sqrt((self.te + self.ti) / mD)

        # What follows has been copy/pasted out of the input file for the diiid case. This is what is used to map from
        # Gkeyll to machine coordinates.
        self.Rdim = 1.7
        self.Zdim = 3.2
        self.Z_axis = 0.013055028
        self.R_axisTrue = 1.6486461
        self.R_axis = 1.6
        self.B_axis = 2 * self.R_axisTrue / self.R_axis
        self.R_LCFSmid = 2.17
        self.Rmid_min = self.R_LCFSmid - 0.1
        self.Rmid_max = self.R_LCFSmid + 0.05
        self.R0 = 0.5 * (self.Rmid_min + self.Rmid_max)
        self.a_mid = self.R_LCFSmid - self.R_axis
        self.minor_r0 = self.R0 - self.R_axis
        self.B0 = self.B_axis * (self.R_axis / self.R0)
        self.qSep = 5.22
        self.sSep = 1.27976219
        self.kappa = 1.35
        self.delta = 0.4
        self.Lz = 2 * np.pi + 1e-8
        self.Lx = self.Rmid_max - self.Rmid_min
        self.xMin, self.xMax = 0., self.Lx
        self.rMin, self.rMax = self.Rmid_min - self.R_axis, self.Rmid_max - self.R_axis
        self.q0 = self.qprofile(self.r_x(0.5 * (self.xMin + self.xMax)))

        # Get the machine coordinates for every point. This takes a long time!
        try:
            XX = np.load(root + "XX.npy")
            YY = np.load(root + "YY.npy")
            ZZ = np.load(root + "ZZ.npy")
            RR = np.load(root + "RR.npy")
            PP = np.load(root + "PP.npy")
            BB = np.load(root + "BB.npy")

        except FileNotFoundError:
            print("Calculating machine coordinates...")
            xx, yy, zz = np.meshgrid(xGrid, yGrid, zGrid)
            XX = np.zeros(xx.shape)
            YY = np.zeros(xx.shape)
            ZZ = np.zeros(xx.shape)
            RR = np.zeros(xx.shape)
            PP = np.zeros(xx.shape)
            BB = np.zeros(xx.shape)
            for i in tqdm(range(0, xx.shape[0])):
                for j in range(0, xx.shape[1]):
                    for k in range(0, xx.shape[2]):
                        coords = self.mapc2p(xx[i, j, k], yy[i, j, k], zz[i, j, k])
                        XX[i, j, k] = coords[0]
                        YY[i, j, k] = coords[1]
                        ZZ[i, j, k] = coords[2]
                        RR[i, j, k] = coords[3]
                        PP[i, j, k] = coords[4]
                        BB[i, j, k] = coords[5]
            np.save(root + "XX.npy", XX)
            np.save(root + "YY.npy", YY)
            np.save(root + "ZZ.npy", ZZ)
            np.save(root + "RR.npy", RR)
            np.save(root + "PP.npy", PP)
            np.save(root + "BB.npy", BB)
        self.XX = XX
        self.YY = YY
        self.ZZ = ZZ
        self.RR = RR
        self.PP = PP
        self.BB = BB

        print("Calculating the magnetic field gradient...")
        gradB = np.gradient(BB, xGrid, yGrid, zGrid)
        self.gradBx = gradB[0]  # Units of T / m
        self.gradBy = gradB[1]  # ''
        self.gradBz = gradB[2]  # Units of T / radian

    def r_x(self, x):
        """
        Go from simulation coordinate to minor radius.
        """
        return x + self.a_mid

    def R(self, r, theta):
        return self.R_axis + r * np.cos(theta + np.arcsin(self.delta) * np.sin(theta))

    def Z(self, r, theta):
        return self.Z_axis + self.kappa * r * np.sin(theta)

    def Bphi(self, R):
        return self.B0 * self.R0 / R

    def dZdr(self, r, theta):
        return self.kappa * np.sin(theta)

    def dZdtheta(self, r, theta):
        return self.kappa * r * np.cos(theta)

    def dRdr(self, r, theta):
        return np.cos(np.arcsin(self.delta) * np.sin(theta) + theta)

    def dRdtheta(self, r, theta):
        return -r * np.sin(np.arcsin(self.delta) * np.sin(theta) + theta) * (np.arcsin(self.delta) * np.cos(theta) + 1)

    def Jr(self, r, theta):

        # df is differentiation of the function. So df(R,1) is differentiation of R wrt to its first argument (r) and
        # df(R,2) is differentiation wrt to its second argument (theta). We just implement functions where I have calculated
        # the derivatives (e.g., WolframAlpha) and not rely on automatic differentiation.
        # return R(r,theta)*(df(R,1)(r,theta)*df(Z,2)(r,theta) - df(Z,1)(r,theta)*df(R,2)(r,theta))
        # print("dRdr = {}".format(dRdr(r, theta)))
        # print("dZdt = {}".format(dZdtheta(r, theta)))
        # print("dZdr = {}".format(dZdr(r, theta)))
        # print("dRdt = {}".format(dRdtheta(r, theta)))
        return self.R(r, theta) * (
                    self.dRdr(r, theta) * self.dZdtheta(r, theta) - self.dZdr(r, theta) * self.dRdtheta(r, theta))

    def dPsidr(self, r, theta):
        def integrand(t): return self.Jr(r, t) / self.R(r, t) ** 2

        # print(r)
        # print(theta)
        integral, _ = integrate.quad(integrand, 0, 2 * np.pi, epsrel=1e-10)
        # print(integral)
        # print(_)
        return self.B0 * self.R_axis / (2 * np.pi * self.qprofile(r)) * integral

    # This one required some rewriting.
    def alpha(self, r, theta, phi):
        def integrand(t):
            return self.Jr(r, t) / self.R(r, t) ** 2

        t = theta
        while t < -np.pi:
            t = t + 2 * np.pi
        while np.pi < t:
            t = t - 2 * np.pi
        if 0 < t:
            # integral, _ =  quad.dblexp(integrand, 0, t, 1e-10)/dPsidr(r,theta)
            integral, _ = integrate.quad(integrand, 0, t, epsabs=1e-10)
            integral = integral / self.dPsidr(r, theta)
        else:
            # integral, _ = -quad.dblexp(integrand, t, 0, 1e-10)/dPsidr(r,theta)
            integral, _ = integrate.quad(integrand, t, 0, epsabs=1e-10)
            integral = integral / self.dPsidr(r, theta)

        return phi - self.B0 * self.R_axis * integral

    def gradr(self, r, theta):
        # return R(r,theta)/Jr(r,theta)*np.sqrt(df(R,2)(r,theta)^2 + df(Z,2)(r,theta)^2)
        return self.R(r, theta) / self.Jr(r, theta) * np.sqrt(
            self.dRdtheta(r, theta) ** 2 + self.dZdtheta(r, theta) ** 2)

    def mapc2p(self, x, y, z):
        r = self.r_x(x)
        R1 = self.R(r, z)
        Z1 = self.Z(r, z)
        phi = -self.q0 / self.minor_r0 * y - self.alpha(r, z, 0)
        X1 = R1 * np.cos(phi)
        Y1 = R1 * np.sin(phi)
        B = self.bmag(x, y, z)
        return X1, Y1, Z1, R1, phi, B

    def qprofile(self, r):
        a = [49.46395467479657, -260.79513158768754, 458.42618139184754, -267.63441353752336]
        return a[0] * (r + self.R_axis) ** 3 + a[1] * (r + self.R_axis) ** 2 + a[2] * (r + self.R_axis) + a[3]

    def bcShiftFunc(self, t, x, y, z):
        """
        This function is used to go from the X simulation coordinate into the length of field line length. In the Gkeyll
        input file this is used to set the boundary condition in the Z direction for the SOL region.
        """
        r = self.r_x(x)
        return self.minor_r0 / self.q0 * self.qprofile(r) * self.Lz

    def bmag(self, x, y, z):
        r = self.r_x(x)
        Bt = self.Bphi(self.R(r, z))
        Bp = self.dPsidr(r, z) / self.R(r, z) * self.gradr(r, z)
        return np.sqrt(Bt**2 + Bp**2)

class Plasma3:

    def __init__(self, root, fstart, fend):
        """
        root (str): Path to the directory containing the background data.
        fstart (int): Starting frame number.
        fend (end): Ending frame number.
        """

        # Load in the plasma background data. This is taken from Tess's nstx-diagnostic.py script.
        step = 1
        Nt = fend - fstart
        Nt //= step
        dt = 0.5e-6

        fp = "{:}d3d-167196-v2".format(root)
        # fp = "{:}d3d-167196".format(root)
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
            d = pg.GData("%s_ion_Temp_%d.bp" % (fp, t))
            dg = pg.GInterpModal(d, 1, 'ms')
            X, tIon = dg.interpolate(0)
            tIon = tIon[:, :, :, 0] / eV
            tIon_tot.append(tIon)

            # Phi eq profile calcuation
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
        phi_tot = np.array(phi_tot)
        elecx_tot = np.array(elecx_tot)
        elecy_tot = np.array(elecy_tot)
        elecz_tot = np.array(elecz_tot)
        self.ne = nElc_tot
        self.te = tElc_tot
        self.ti = tIon_tot
        self.vp = phi_tot
        self.er = elecx_tot
        self.ep = elecy_tot
        self.ez = elecz_tot
        self.r = x
        self.p = y
        self.z = z
        self.dt = dt
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

# Initialize background. Choose from which case here.
case = "diiid-167196"
if case == "nstx":
    plasma = Plasma("/Users/zamperini/Documents/d3d_work/gkyl_files/nstx-neut/", 200, 450)
elif case == "diiid":
    fstart = 500
    fend = 654
    plasma = Plasma2(tstart=500, tend=654)
    plasma.r = plasma.r - 0.10   # Convert to R-Rsep
elif case == "diiid-167196":
    fstart = 600
    fend = 999
    plasma = Plasma3("/Users/zamperini/Documents/d3d_work/gkyl_files/d3d-167196-v2/", 600, 999)  # 0-999
else:
    print("Error! case must be one of nstx or diiid.")
    plasma = None
    sys.exit()

# Option to modify the background by interpolating additional frames in between each.
interp_frames = True
if interp_frames:
    nadd = 5  # Number of frames to add between each.
    nframes = plasma.ne.shape[0]
    old_dt = plasma.dt
    new_dt = old_dt / (nadd + 1)
    print("Interpolating frames...")
    arrs = [plasma.ne, plasma.te, plasma.ti, plasma.er, plasma.ep, plasma.ez, plasma.cs]
    tmp_arrs = []
    for i in tqdm(range(0, len(arrs))):
        # print("Array {}/{}".format(i+1, len(arrs)))
        tmp_arr = np.zeros((((nframes - 1) * nadd + 1), arrs[i].shape[1], arrs[i].shape[2], arrs[i].shape[3]))
        for f in range(0, nframes-1):
            y0 = arrs[i][f]
            y1 = arrs[i][f + 1]
            m = (y1 - y0) / old_dt
            tmp_arr[f * nadd] = y0
            tmp_arr[(f+1) * nadd] = y1
            for g in range(1, nadd+1):
                # arrs[i] = np.insert(arrs[i], (f * nadd) + g, m * new_dt * g + y0, axis=0)
                tmp_arr[f * nadd + g] = m * new_dt * g + y0
        tmp_arrs.append(tmp_arr)
    plasma.ne = tmp_arrs[0]
    plasma.te = tmp_arrs[1]
    plasma.ti = tmp_arrs[2]
    plasma.er = tmp_arrs[3]
    plasma.ep = tmp_arrs[4]
    plasma.ez = tmp_arrs[5]
    plasma.cs = tmp_arrs[6]
    plasma.dt = new_dt

# Bounds.
rmin = plasma.r.min()
rmax = plasma.r.max()
pmin = plasma.p.min()
pmax = plasma.p.max()
midz = len(plasma.z) // 2

# A synthetic (and real) collector probe. Since we are in the R, P plane, we are "looking at" the faces of the CP. We
# ignore parallel transport, so technically there is no way for the ions to reach the probe faces; they would all hit
# the edges. This is on shaky ground.
include_cp = False
cp_a2 = pd.read_excel("/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx")
cp_r = np.linspace(cp_a2["R D (cm)"].min(), cp_a2["R D (cm)"].max(), 100) / 100
cp_deps1 = np.zeros(cp_r.shape)
cp_deps2 = np.zeros(cp_r.shape)

# Impurity simulation settings.
imp_mode = "montecarlo"
sol_loss = False
mem_lite = False  # Do not count particle history to avoid huge memory load.
imp_Z = 15
imp_mass = 184.84 * 1.66E-27  # amu x kg
imps = []
imp_dens = np.zeros((plasma.ne.shape[1], plasma.ne.shape[2], plasma.ne.shape[3]))
imp_vrs = np.zeros((plasma.ne.shape[1], plasma.ne.shape[2], plasma.ne.shape[3]))
imp_cnts = np.zeros((plasma.ne.shape[1], plasma.ne.shape[2], plasma.ne.shape[3]))
print("Impurity starting location: {}".format(imp_mode))
if imp_mode == "grid":

    # Initialize equally spaced impurities in a 2D grid.
    nrows = 25
    ncols = 25
    # rstarts = np.linspace(rmin, rmax, ncols)
    if case == "nstx":
        rstarts = np.linspace(1.30, 1.34, ncols)
    elif case == "diiid":
        rstarts = np.linspace(-0.01, 0.01, ncols)
        # rstarts = np.linspace(rmin, rmax, ncols)
    elif case == "diiid-167196":
        rstarts = np.linspace(rmin, rmax, ncols)
    else:
        rstarts = None
    pstarts = np.linspace(pmin, pmax, nrows)
    for row in range(0, nrows):
        for col in range(0, ncols):
            # Create impurity ions. Time step must equal plasma.dt!!! That idea didn't pan out yet.
            imp = Impurity(imp_Z, imp_mass, rstarts[col], pstarts[row], 0.0, dt=plasma.dt)
            imps.append(imp)

elif imp_mode == "montecarlo":
    nimps = 10000
    if case == "nstx":
        inj_rmin = 1.31
        inj_rmax = 1.34
    elif case == "diiid":
        #inj_rmin = 0.0  # Separatrix
        #inj_rmax = 0.0
        inj_rmin = 0.01  # Couple cm outside separatrix
        inj_rmax = 0.01
        # inj_rmin = rmin
        # inj_rmax = rmax
    elif case == "diiid-167196":
        # inj_rmin = rmin
        # inj_rmax = rmax
        inj_rmin = 2.32
        inj_rmax = 2.32
        # inj_rmin = 2.38
        # inj_rmax = 2.38
    inj_rs = inj_rmin + (inj_rmax - inj_rmin) * np.random.random(nimps)
    inj_ps = pmin + (pmax - pmin) * np.random.random(nimps)

    # To add soem noise to the simulation video, don't start all the impurities at the same time. Allow them to come in
    # at random times throughout the simulation.
    inj_fs = np.array(fstart + (fstart - fend) * np.random.random(nimps), dtype=int)
    # inj_fs = np.full(inj_rs.shape, fstart)  # Start all at the beginning like usual.

    # inj_zs = 2 * np.pi * np.random.random(nimps)
    inj_zs = np.full(nimps, plasma.z[midz])
    for i in range(0, nimps):
        imp = Impurity(imp_Z, imp_mass, inj_rs[i], inj_ps[i], inj_zs[i], dt=plasma.dt, fbirth=inj_fs[i])
        imps.append(imp)

# Go through one timeframe at a time, updating the impurities positions based on just the ExB force.
bt = 2.0  # Assuming constant B field for now.
bp = 0.5
b = 2.0
e = 1.609e-19
omega_c = imp_Z * e * b / imp_mass
if case == "nstx":
    rsep = 1.32
    lsol = 100.0
elif case == "diiid":
    rsep = 0.0
    lsol = 100.0
elif case == "diiid-167196":
    rsep = 2.30
    lsol = 100.0
else:
    print("Error? Case name not recognized...")
    rsep = 0.0
    lsol = 0.0

prev_exb_drft_r = 0.0
prev_exb_drft_p = 0.0
prev_exb_drft_z = 0.0
prev_p = 0.0
print("Following impurities...")
for f in tqdm(range(0, plasma.ne.shape[0])):

    # Index the electric field in the (R,P) plane.
    er = plasma.er[f, :, :, :]
    ep = plasma.ep[f, :, :, :]
    ez = plasma.ez[f, :, :, :]
    cs = plasma.cs[f, :, :, :]
    ti = plasma.ti[f, :, :, :]

    if f == 0:
        tstart = 0
    else:
        tstart = (f - 1) * plasma.dt
    tend = f * plasma.dt

    for imp in imps:

        # If ion not born yet, skip (by putting in nans).
        if (imp.fbirth - fstart) > f:
            if not mem_lite:
                imp.rhist.append(np.nan)
                imp.phist.append(np.nan)
                imp.vrhist.append(np.nan)
                imp.vphist.append(np.nan)
                continue

        # Calculate number of steps we are going to take.
        # nsteps = int((tend - tstart) / imp.dt)
        nsteps = 1
        t_in_sol = 0.0
        skip_pol_drft = True
        for step in range(0, nsteps):

            # Killed options below.
            if imp.killed:

                # 8/15/23: Instead of just killing the impurity, reset it randomly. This has only been thought about in
                # a Monte Carlo sense. We skip the polarization drift at first since it is starting again from rest, and
                # we don't want to give a huge value for the time derivative.
                if imp_mode == "montecarlo":
                    imp.r = inj_rmin + (inj_rmax - inj_rmin) * np.random.random()
                    imp.p = pmin + (pmax - pmin) * np.random.random()
                    imp.vr = 0.0
                    imp.vp = 0.0
                    skip_pol_drft = True
                    imp.killed = False

                # Otherwise don't do anything if this impurity was killed.
                else:
                    if not mem_lite:
                        imp.rhist.append(np.nan)
                        imp.phist.append(np.nan)
                        imp.vrhist.append(np.nan)
                        imp.vphist.append(np.nan)
                        imp.ephist.append(np.nan)
                    continue

            # Index the cell we're in.
            ridx = np.argmin(np.abs(imp.r - plasma.r))
            pidx = np.argmin(np.abs(imp.p - plasma.p))
            zidx = np.argmin(np.abs(imp.z - plasma.z))
            imp_dens[ridx, pidx, zidx] += imp.dt
            imp_cnts[ridx, pidx, zidx] += 1

            # Magnetic moment, mu = m * v_perp^2 / 2B and v_perp^2 ~ 2Ti / m, thus mu = m * (2Ti / m) / 2B = Ti / B.
            vperp = np.sqrt(2 * ti[ridx, pidx, zidx] * e / imp_mass)  # Units of m/s
            mag_mom = ti[ridx, pidx, zidx] / plasma.BB[ridx, pidx, zidx] * e  # Units of J / T.

            # ExB drift for each direction. In the z direction the units are actually rad/s (since when we took the
            # gradient of E in the z direction it used the units of z, radians, as the unit).
            exb_drft_r = ep[ridx, pidx, zidx] / plasma.BB[ridx, pidx, zidx]
            exb_drft_p = -er[ridx, pidx, zidx] / plasma.BB[ridx, pidx, zidx]
            exb_drft_z = 0.0

            # The curvature drift, where vpar ~ vthermal. This is in the binormal (p) direction. For the
            # diiid case the values are not the actual R values, so I index it from the derived values. The parallel
            # velocity is approximated as the sound speed, otherwise we have to solve a differential equation that
            # requires the gradient of B along the field line, dB/dl, which is not readily available.
            if case == "diiid":
                curve_drft_p = -cs[ridx, pidx, zidx] ** 2 / omega_c / plasma.RR[ridx, pidx, zidx]
            else:
                curve_drft_p = -cs[ridx, pidx, zidx] ** 2 / omega_c / imp.r

            # The polarization drift for each direction. Uses previous ExB drift value to calculate the needed
            # derivative.
            if f == 0 or skip_pol_drft:
                pol_drft_r = 0.0
                pol_drft_p = 0.0
                pol_drft_z = 0.0
                skip_pol_drft = False
            else:
                dexb_dt_r = (exb_drft_r - prev_exb_drft_r) / imp.dt
                dexb_dt_p = (exb_drft_p - prev_exb_drft_p) / imp.dt
                # dexb_dt_z = (exb_drft_z - prev_exb_drft_z) / imp.dt

                # Note: The r and p here are intentional because the polarization drift is b x d(vExB)/dt, where b is
                # the direction of the magnetic field (not the whole vector, just b hat).
                pol_drft_r = -prev_exb_drft_p / omega_c
                pol_drft_p = prev_exb_drft_r / omega_c
                pol_drft_z = 0.0

            # The grad-B drift. This is the same whether we say gradBz = 0 or not; the gradBz terms are all with Bx and
            # By terms, which zero those terms out anyways.
            gradb_drft_r = -vperp ** 2 / omega_c * plasma.gradBy[ridx, pidx, zidx] / plasma.BB[ridx, pidx, zidx]
            gradb_drft_p = vperp ** 2 / omega_c * plasma.gradBx[ridx, pidx, zidx] / plasma.BB[ridx, pidx, zidx]

            # Assign the particle velocity to the local drifts.
            imp.vr = exb_drft_r + pol_drft_r + gradb_drft_r
            imp.vp = exb_drft_p + pol_drft_p + gradb_drft_p + curve_drft_p
            imp.vz = 0.0
            # imp.vr = exb_drft_r
            # imp.vp = exb_drft_p
            # imp.vz = exb_drft_z
            imp.r = imp.r + imp.vr * imp.dt
            imp.p = imp.p + imp.vp * imp.dt
            imp.z = imp.z + imp.vz * imp.dt

            # If we are out of the simulation bound we need to kill the particle.
            if imp.r > rmax or imp.r < rmin:
                # print("Killed: {:.3f}  {:.3f}".format(imp.r, imp.p))
                imp.killed = True

                # If mem_lite is on, save this final position since it is needed below. Minus to not include the
                # contribution that was just added.
                if mem_lite:
                    imp.rhist.append(imp.r - imp.vr * imp.dt)
                    imp.phist.append(imp.p - imp.vp * imp.dt)

            # Condition if we are in the SOL and reached the edge of the field line (where a hypothetical limiter is).
            elif case == "diiid" and imp.r > rsep and (imp.z < 0.0 or imp.z > 2 * np.pi):
                imp.killed = True

                # If mem_lite is on, save this final position since it is needed below. Minus to not include the
                # contribution that was just added.
                if mem_lite:
                    imp.rhist.append(imp.r - imp.vr * imp.dt)
                    imp.phist.append(imp.p - imp.vp * imp.dt)

            else:
                imp.killed = False

            # Assume repeating bounds in the poloidal direction, so once a particle passes out of bounds just assign
            # it as wrapping around on the other poloidal side.
            if imp.p > pmax:
                imp.p = pmin + (imp.p - pmax)
            if imp.p < pmin:
                imp.p = pmax + (imp.p - pmin)

            # Add to the particle history lists. These are the values at the end of the steps.
            if not mem_lite:
                imp.rhist.append(imp.r)
                imp.phist.append(imp.p)
                imp.zhist.append(imp.z)
                imp.vrhist.append(imp.vr)
                imp.vphist.append(imp.vp)
                imp.vzhist.append(imp.vz)
                imp.ephist.append(ep[ridx, pidx, zidx])

            # If particle is in the SOL, randomly see if we lose that particle using L/cs for the characteristic loss
            # time. Separatrix values are hardcoded.
            if sol_loss:
                if imp.r > rsep:
                    t_in_sol += imp.dt
                else:
                    t_in_sol = 0.0
                loss_time = lsol / cs[ridx, pidx]

                # Use rejection method with exponential for if the particle is killed off.
                if t_in_sol > np.exp(-t_in_sol / loss_time):
                    imp.killed = True
                    if mem_lite:
                        imp.rhist.append(imp.r - imp.vr * imp.dt)
                        imp.phist.append(imp.p - imp.vp * imp.dt)

            # See what deposition on a synthetic CP could look like, assuming it as at P = 0.
            if include_cp:
                if imp.r > cp_r.min():

                    # Crossed from "bottom"
                    if prev_p < 0 and imp.p > 0:
                        cp_idx = np.argmin(np.abs(imp.r - cp_r))
                        cp_deps1[cp_idx] += 1
                        imp.killed = True
                        if mem_lite:
                            imp.rhist.append(imp.r - imp.vr * imp.dt)
                            imp.phist.append(imp.p - imp.vp * imp.dt)

                    # Crossed from "top"
                    if prev_p > 0 and imp.p < 0:
                        cp_idx = np.argmin(np.abs(imp.r - cp_r))
                        cp_deps2[cp_idx] += 1
                        imp.killed = True
                        if mem_lite:
                            imp.rhist.append(imp.r - imp.vr * imp.dt)
                            imp.phist.append(imp.p - imp.vp * imp.dt)

            # Add to array, so we can calculate an average later.
            imp_vrs[ridx, pidx, zidx] += imp.vr

            # Save for the derivative calculation.
            prev_exb_drft_r = exb_drft_r
            prev_exb_drft_p = exb_drft_p
            prev_exb_drft_z = exb_drft_z
            prev_p = imp.p

# Plot comparing how far a particle traveled to the amount of time spent in negative Ep zones.
pos_min_neg = []
dist_trav = []
for imp in imps:
    nonan = ~np.isnan(imp.ephist)
    neg_ep_time = np.sum(np.array(imp.ephist)[nonan] < 0.0) * imp.dt
    pos_ep_time = np.sum(np.array(imp.ephist)[nonan] > 0.0) * imp.dt
    dist_trav.append(np.array(imp.rhist)[nonan][-1] - np.array(imp.rhist)[nonan][0])
    pos_min_neg.append(pos_ep_time - neg_ep_time)
pos_min_neg = np.array(pos_min_neg) * 1e6
dist_trav = np.array(dist_trav) * 100
fig, ax = plt.subplots(figsize=(7, 6))
ax.axvline(0.0, color="k", zorder=5)
ax.axhline(0.0, linestyle="--", color="k", zorder=5)
ax.scatter(pos_min_neg, dist_trav, color="tab:red", edgecolors="k", zorder=15)
ax.set_xlabel(r"$\mathdefault{\Delta t\ (\mu s)}$", fontsize=16)
ax.set_ylabel(r"$\mathdefault{\Delta R\ (cm)}$", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.show()


# Calculate average radial velocity in each cell.
for i in range(0, imp_vrs.shape[0]):
    for j in range(0, imp_vrs.shape[1]):
        for k in range(0, imp_vrs.shape[2]):
            if imp_cnts[i, j, k] == 0:
                imp_vrs[i, j, k] = 0.0
            else:
                imp_vrs[i, j, k] = imp_vrs[i, j, k] / imp_cnts[i, j, k]

# Normalize and divide by area of box (and 1 m-tor) to get m-3.
imp_dens = imp_dens / imp_dens.sum()
for i in range(0, imp_dens.shape[0]):
    for j in range(0, imp_dens.shape[1]):
        for k in range(0, imp_dens.shape[2]):

            # dr and dp are in units, fine.
            if i == imp_dens.shape[0] - 1:
                dr = plasma.r[i] - plasma.r[i - 1]
            else:
                dr = plasma.r[i + 1] - plasma.r[i]
            if j == imp_dens.shape[1] - 1:
                dp = plasma.p[j] - plasma.p[j - 1]
            else:
                dp = plasma.p[j + 1] - plasma.p[j]

            # dz is naturally in radians, we want it in meters. The field line length is approximated as L = pi * R * q
            # (Eq. 1.7 in Peter's book, length = 2L). For now use approximation for a circular cross-section such that
            # every dz = constant.
            if case == "diiid":
                L = np.pi * plasma.RR[i, j, k] * plasma.qprofile(plasma.r_x(plasma.r[i]))
                # if k == imp_dens.shape[2] - 1:
                #     dz = (plasma.z[k] - plasma.z[k - 1]) / (2 * np.pi) * L
                # else:
                #     dz = (plasma.z[k + 1] - plasma.z[k]) / (2 * np.pi) * L
                dz = 1.0
            elif case in ["nstx", "diiid-167196"]:
                dz = plasma.z[1] - plasma.z[0]

            imp_dens[i, j, k] = imp_dens[i, j, k] / (dr * dp * dz)

imp_dens[imp_dens == 0.0] = imp_dens[imp_dens != 0.0].min()

if case == "diiid":
    fig, ax = plt.subplots()
    # ax.tricontourf(imp_dens_r.flatten(), imp_dens_z.flatten(), imp_cnts.flatten(), locator=ticker.LogLocator())
    ax.tricontourf(plasma.RR.flatten(), plasma.ZZ.flatten(), imp_cnts.flatten())
    ax.set_aspect("equal")
    ax.set_xlim(0.75, 2.5)
    ax.set_ylim(-1.1, 1.1)
    fig.tight_layout()
    fig.show()

# For each radial bin (e.g., every cm) calculate the average distance traveled by ions that started there.
bins = np.linspace(rmin, rmax, 15)
bin_sums = np.zeros(bins.shape)
bin_counts = np.zeros(bins.shape)
bin_vals = {b: [] for b in bins}
for imp in imps:
    for i in range(0, len(bins) - 1):
        if bins[i] <= imp.rhist[0] < bins[i + 1]:
            # Calculate radial distance traveled. Add to bin sum.
            raddist = np.array(imp.rhist)[~np.isnan(imp.rhist)][-1] - imp.rhist[0]
            bin_sums[i] += raddist
            bin_counts[i] += 1
            bin_vals[bins[i]].append(raddist)
# avg_raddist = bin_sums / bin_counts
avg_raddist = np.array([bin_sums[i] / bin_counts[i] if bin_counts[i] > 0 else np.nan for i in range(0, len(bin_counts))])
med_raddist = np.array([np.median(bin_vals[b]) for b in bins])

# Pull out a radial "density" profile by averaging along the binormal direction.
# rad_imp_dens = imp_dens[:, imp_dens.shape[1] // 2]
imp_dens[imp_dens == 0] = np.nan
rad_imp_dens = savgol_filter(np.nanmean(imp_dens[:, :, midz], axis=1), 11, 3)

# Extract a comparable diffusion coefficient from the fluxes.
rad_imp_vrs = np.nanmean(imp_vrs[:, :, midz], axis=1)
rad_imp_gmma = rad_imp_vrs * rad_imp_dens
rad_imp_dndr = np.zeros(rad_imp_dens.shape)
delth = plasma.r[1] - plasma.r[0]
for i in range(1, len(rad_imp_dndr)):
    if i < len(rad_imp_dndr) - 1:
        rad_imp_dndr[i] = (rad_imp_dens[i + 1] - rad_imp_dens[i]) / (plasma.r[i + 1] - plasma.r[i])
        # rad_imp_dndr[i] = (rad_imp_dens[i + 1] - rad_imp_dens[i - 1]) / (2 * delth)
rad_imp_diff = np.zeros(rad_imp_gmma.shape)
for i in range(0, len(rad_imp_diff)):
    if rad_imp_dndr[i] != 0:
        rad_imp_diff[i] = -rad_imp_gmma[i] / rad_imp_dndr[i]

segs = []
colors = []
i = 0
data_type = "velocity"  # "velocity", "gamma" or "density"
if data_type == "gamma":
    data = rad_imp_gmma
elif data_type == "density":
    data = rad_imp_dens
elif data_type == "velocity":
    data = rad_imp_vrs
else:
    print("Error in data_type: {}".format(data_type))
    data = None
while i < len(data):
    tmp_seg = []
    tmp_r = []
    sign = np.sign(data[i])
    val_sign = np.sign(data[i])
    while val_sign == sign:
        tmp_seg.append(data[i] * val_sign)
        tmp_r.append(plasma.r[i])
        i += 1
        if i == len(data): break
        val_sign = np.sign(data[i])

    # Add the next value so the lines connect.
    if i < len(data) - 1:
        tmp_r.append(plasma.r[i + 1])
        tmp_seg.append(data[i + 1] * val_sign)

    segs.append(list(zip(tmp_r, tmp_seg)))
    if sign == -1:
        colors.append("tab:purple")
    else:
        colors.append("tab:red")

split_velocities = False
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
if case == "nstx":
    ax1.set_xlim([1.27, 1.42])
    ax1.text(1.325, 1.0, "Quasi-\nseparatrix", fontsize=12, horizontalalignment="center")
elif case == "diiid":
    # ax1.text(0.10, 1.0, "Separatrix", fontsize=12, horizontalalignment="center")
    if data_type == "gamma":
        ax1.set_ylim([1e1, 1e5])
    elif data_type == "density":
        ax1.set_ylim([5e0, 1e2])
    ax1.axvline(0.0, color="k", linestyle="--", lw=3)
    ax2.axvline(0.0, color="k", linestyle="--", lw=3)
    pass
elif case == "diiid-167196":
    ax1.set_xlim([rmin, rmax])
    ax1.axvline(2.31, color="k", linestyle="--", lw=3)
    ax2.axvline(2.31, color="k", linestyle="--", lw=3)
    # ax1.text(2.33, 225, "Quasi-\nseparatrix", fontsize=12, horizontalalignment="center")
if imp_mode == "montecarlo":
    # ax1.axvspan(inj_rmin, inj_rmax, color="grey", alpha=0.3)
    # ax2.axvspan(inj_rmin, inj_rmax, color="grey", alpha=0.3)
    pass
if split_velocities:
    line_segs1 = LineCollection(segs, linewidths=3, colors="k", linestyle="solid")
    line_segs2 = LineCollection(segs, linewidths=2, colors=colors, linestyle="solid")
    ax1.add_collection(line_segs1)
    ax1.add_collection(line_segs2)
else:
    ax1.axhline(0.0, color="k", lw=2)
    ax1.plot(plasma.r, data, color="k", lw=3)
    ax1.plot(plasma.r, data, color="tab:red", lw=2)
ax1.set_xlabel("Radial (m)", fontsize=12)
if data_type == "gamma":
    ax1.set_ylabel("Impurity Flux (arbitrary)", fontsize=12)
    ax1.set_yscale("log")
elif data_type == "denisty":
    ax1.set_ylabel("Impurity Density (arbitrary)", fontsize=12)
    ax1.set_yscale("log")
else:
    ax1.set_ylabel("Radial Velocity (m/s)", fontsize=12)
    ax1.set_ylim([-100, 300])
ax1.grid(zorder=-1)

ax2.plot(plasma.r, rad_imp_diff, color="k", lw=3)
ax2.plot(plasma.r, rad_imp_diff, color="tab:red", lw=2)
# ax2.scatter(plasma.r, rad_imp_diff, color="tab:red", edgecolor="k", zorder=50)
ax2.set_ylabel(r"$\mathdefault{D}^{eff}_{\perp}\ \mathdefault{(m^2/s)}$", fontsize=12)
ax2.set_xlabel("Radial (m)", fontsize=12)
ax2.grid(zorder=-1)
ax2.set_yscale("log")
if case == "nstx":
    ax2.set_ylim([0.01, 30])
elif case == "diiid":
    ax2.set_ylim(0.001, 1000)
elif case == "diiid-167196":
    ax2.set_ylim([0.1, 100])
fig.tight_layout()
fig.show()

# Plot considering a synthetic probe measuring blobs.
sig_locr = np.linspace(2.32, 2.38, 5)
sig_locp = np.zeros(len(sig_locr))
sig_t = np.arange(0, plasma.ne.shape[0] * imp.dt, imp.dt)
sig_ne = np.zeros((sig_locr.shape[0], plasma.ne.shape[0]))
sig_ep = np.zeros((sig_locr.shape[0], plasma.ne.shape[0]))
for f in tqdm(range(0, plasma.ne.shape[0])):
    ne = plasma.ne[f]
    ep = plasma.ep[f]
    for i in range(0, len(sig_locr)):
    #for r, p in zip(sig_locr, sig_locp):
        r = sig_locr[i]
        p = sig_locp[i]
        idx_r = np.argmin(np.abs(r - plasma.r))
        idx_p = np.argmin(np.abs(p - plasma.p))
        sig_ne[i, f] = ne[idx_r, idx_p, midz]
        sig_ep[i, f] = ne[idx_r, idx_p, midz]
rms_ne = np.array([np.sqrt(np.mean(np.square(sig_ne[i]))) for i in range(0, sig_ne.shape[0])])

colors = ["C{}".format(i) for i in range(0, 10)]
fig, ax1 = plt.subplots()
for i in range(0, sig_ne.shape[0]):
    ax1.plot(sig_t * 1e6, sig_ne[i], label="{:.3f}".format(sig_locr[i]), color=colors[i])
    ax1.axhline(rms_ne[i], linestyle="--", color=colors[i])
ax1.set_xlabel("Time (us)")
ax1.legend()
fig.tight_layout()
fig.show()

# Plot of a synthetic collector probe. See comment above, this is all on shaky grounds.
if include_cp:
    cp_a2_itf = cp_a2["W Areal Density D (1e15 W/cm2)"] / cp_a2["W Areal Density D (1e15 W/cm2)"].max()
    cp_a2_otf = cp_a2["W Areal Density U (1e15 W/cm2)"] / cp_a2["W Areal Density U (1e15 W/cm2)"].max()
    cp_deps1_norm = cp_deps1 / cp_deps1.max()
    cp_deps2_norm = cp_deps1 / cp_deps2.max()
    cp_min = cp_a2_itf[cp_a2_itf > 0].min()
    fig, ax1 = plt.subplots()
    ax1.plot(cp_r, cp_deps1_norm)
    ax1.plot(cp_r, cp_deps2_norm)
    ax1.scatter(cp_a2["R D (cm)"] / 100, cp_a2_itf)
    ax1.scatter(cp_a2["R U (cm)"] / 100, cp_a2_otf)
    ax1.set_yscale("log")
    fig.tight_layout()
    fig.show()

# Some numbers for the plots below.
ep_min = plasma.ep[:, :, :, midz].min()
ep_max = plasma.ep[:, :, :, midz].max()
ep_lim = np.max(np.abs([ep_min, ep_max])) / 2
plot_ep = np.clip(plasma.ep, -ep_lim, ep_lim)
time_label = 0.0
scat_cmap = matplotlib.cm.get_cmap("PiYG")
max_raddist = 0.03

# Plot summary. First assemble arrays of the particle position histories.
back_data = "ep"
fontsize = 12
figsize = (7, 6)
rhists = np.array([imp.rhist for imp in imps])
phists = np.array([imp.phist for imp in imps])
fig, ax = plt.subplots(figsize=figsize)
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
if back_data == "ne":
    cont = ax.contourf(plasma.r, plasma.p, plasma.ne[0, :, :, midz].T, zorder=10)
elif back_data == "ep":
    cont = ax.contourf(plasma.r, plasma.p, plot_ep[0, :, :, midz].T, zorder=10,
                       levels=np.linspace(-ep_lim, ep_lim, 11), cmap="coolwarm")
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{E}_{\perp}$ (V/m)", fontsize=fontsize)

if imp_mode == "grid":
    scat = ax.scatter(rhists[:, 0], phists[:, 0], zorder=20, color="k", edgecolors="k")
elif imp_mode == "montecarlo":
    # npoints = 250
    npoints = 500
    scat = ax.scatter(rhists[:npoints, 0], phists[:npoints, 0], zorder=20, color="k", edgecolors="k")
if case == "nstx":
    ax.set_xlim([1.26, 1.42])
    ax.set_ylim([-0.14, 0.14])
elif case == "diiid":
    ax.axvline(0.0, color="k", zorder=19)
    print("Note: Only showing part of domain!!!")
    ax.set_xlim([-0.03, 0.05])
if case == "diiid":
    ax.set_xlabel("R-Rsep (m)", fontsize=fontsize)
else:
    ax.set_xlabel("Radial (m)", fontsize=fontsize)
ax.set_ylabel("Binormal (m)", fontsize=fontsize)
ax.set_title("{:7.2f} us".format(time_label), fontsize=fontsize)
fig.tight_layout()


def update(frame):
    # Update contour.
    ax.clear()
    if back_data == "ne":
        cont = ax.contourf(plasma.r, plasma.p, plasma.ne[frame, :, :, midz].T, zorder=10)
    elif back_data == "ep":
        cont = ax.contourf(plasma.r, plasma.p, plot_ep[frame, :, :, midz].T, zorder=10,
                           levels=np.linspace(-ep_lim, ep_lim, 11), cmap="coolwarm")

    # Update scatter. Markers are colored by their radial distance from their starting location.
    if imp_mode == "grid":
        raddist = (rhists[:, frame] - rhists[:, 0] + max_raddist) / (max_raddist + max_raddist)
        colors = scat_cmap(raddist)
        scat = ax.scatter(rhists[:, frame], phists[:, frame], zorder=20, color=colors, edgecolors="k")
    elif imp_mode == "montecarlo":
        raddist = (rhists[:npoints, frame] - rhists[:npoints, 0] + max_raddist) / (max_raddist + max_raddist)
        colors = scat_cmap(raddist)
        scat = ax.scatter(rhists[:npoints, frame], phists[:npoints, frame], zorder=20, color=colors, edgecolors="k")

    cax.cla()
    cbar = fig.colorbar(cont, cax=cax, extend="both")
    cbar.set_label(r"$\mathdefault{E}_{\perp}$ (V/m)", fontsize=fontsize)
    if case == "nstx":
        ax.set_xlim([1.26, 1.42])
        ax.set_ylim([-0.14, 0.14])
    elif case == "diiid":
        ax.axvline(0.10, color="k", zorder=19)
        ax.set_xlim([-0.03, 0.05])
    if case == "diiid":
        ax.set_xlabel("R-Rsep (m)", fontsize=fontsize)
    else:
        ax.set_xlabel("Radial (m)", fontsize=fontsize)
    ax.set_ylabel("Binormal (m)", fontsize=fontsize)
    ax.set_title("{:7.2f} us".format(imp.dt * frame * 1e6), fontsize=fontsize)
    fig.tight_layout()
    return scat, cont

print("Creating animation...")
start = time.time()
interval = 80  # Usual value.
# interval = 40
ani = animation.FuncAnimation(fig=fig, func=update, frames=rhists.shape[1] - 1, interval=interval)
end = time.time()
print("Done (took {:.0f} seconds)".format(end - start))
if case == "nstx":
    fname = "nstx-neut-imps"
elif case == "diiid":
    fname = "d3d-posD04-k135-3xv2-imps"
elif case == "diiid-167196":
    fname = "diiid-167196-imps"

start = time.time()
# ani.save(filename="{:}.html".format(fname), writer="html")
end = time.time()
print("HTML saved (took {:.0f} seconds)".format(end - start))
start = time.time()
# ani.save(filename="{:}.gif".format(fname), writer="pillow")
end = time.time()
print("GIF saved (took {:.0f} seconds)".format(end - start))
start = time.time()
if case == "diiid-167196":
    writervideo = animation.FFMpegWriter(fps=25)
elif case == "diiid":
    writervideo = animation.FFMpegWriter(fps=15)
ani.save('{:}.mp4'.format(fname), writer=writervideo)
end = time.time()
print("MP4 saved (took {:.0f} seconds)".format(end - start))
# plt.show()