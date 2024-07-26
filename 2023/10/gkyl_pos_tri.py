import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import sys
from tqdm import tqdm


mD = 931.49e6 / 3e8 ** 2

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


# Load data for indicated time frame.
fstart = 500
fend = 654
plasma = Plasma2(tstart=500, tend=654)
plasma.r = plasma.r - 0.10   # Convert to R-Rsep

plot_idx = 0
fig, ax = plt.subplots()
# ax.tricontourf(imp_dens_r.flatten(), imp_dens_z.flatten(), imp_cnts.flatten(), locator=ticker.LogLocator())
ax.tricontourf(plasma.RR.flatten(), plasma.ZZ.flatten(), plasma.ne[plot_idx].flatten())
ax.set_aspect("equal")
ax.set_xlim(0.75, 2.5)
ax.set_ylim(-1.1, 1.1)
fig.tight_layout()
fig.show()