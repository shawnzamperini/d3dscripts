import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import fsolve


class EngelhardtModel:

    def __init__(self):
        pass

    def set_geometry(self, riz, rlim, rwall):

        self.r0 = -0.05
        self.rsep = 0.0
        self.riz = riz
        self.rlim = rlim
        self.rwall = rwall
        self.rprof = np.linspace(self.r0, self.rwall+0.01, 250)

    def set_ne(self, nesep, lambda_ne):
        neprof = nesep * np.exp(-self.rprof / lambda_ne)

        # Values in the core just cap at the separatrix value.
        neprof[self.rprof <= 0] = nesep

        fne = interp1d(self.rprof, neprof)
        self.fne = fne

        # Assign ne for each region based off innermost value.
        ne = np.array([nesep])  # Core value.
        for r in [self.rsep, self.riz, self.rlim]:
            ne = np.append(ne, fne(r))
        self.ne = ne

    def set_te_ti(self, tesep, lambda_te, timult=3):
        teprof = tesep * np.exp(-self.rprof / lambda_te)
        teprof[self.rprof <= 0] = tesep
        fte = interp1d(self.rprof, teprof)
        fti = interp1d(self.rprof, teprof * timult)
        self.fte = fte
        self.fti = fti
        self.timult = timult

        # Assign Te, Ti for each region based off innermost value.
        te = np.array([tesep]); ti = np.array([tesep * timult])
        for r in [self.rsep, self.riz, self.rlim]:
            te = np.append(te, fte(r))
            ti = np.append(ti, fti(r))
        self.te = te
        self.ti = ti

    def load_mm_model(self):
        # First we want to load in all the yields.
        mmpath = "/Users/zamperini/My Drive/Research/Documents/2022/02/mixed_material.xlsx"
        columns = ["D-Si_E", "D-Si_Y", "D-Si_Ech", "D-Si_Ych25", "D-Si_Ych300",
          "D-Si_Ych600", "D-C_E", "D-C_Y", "D-C_Ech", "D-C_Ychsurf", "D-C_Ychphys",
          "D-SiC,C_E", "D-SiC,C_Y", "D-SiC,C_Ech", "D-SiC,C_Ychsurf", "D-SiC,Si_E",
          "D-SiC,Si_Y", "D-SiC,Si_Ech", "D-SiC,Si_Ychsurf", "C-Si_E",
          "C-Si_Y", "C-C_E", "C-C_Y", "C-SiC,C_E", "C-SiC,C_Y", "C-SiC,Si_E",
          "C-SiC,Si_Y", "Si-C_E", "Si-C_Y", "Si-Si_E", "Si-Si_Y", "Si-SiC,C_E",
          "Si-SiC,C_Y", "Si-SiC,Si_E", "Si-SiC,Si_Y"]
        mm = pd.read_excel(mmpath, skiprows=6, header=None, names=columns,
            usecols="A:F,H:L,N:Q,S:V,X,Y,AA,AB,AD,AE,AG,AH,AJ,AK,AM,AN,AP,AQ,AS,AT")

        # Create interpolation functions for each yield.
        self.Y_D_Si = interp1d(mm["D-Si_E"], mm["D-Si_Y"], fill_value=0, bounds_error=False)
        self.Y_D_Si_ch25 = interp1d(mm["D-Si_Ech"], mm["D-Si_Ych25"], fill_value=0, bounds_error=False)
        self.Y_D_Si_ch300 = interp1d(mm["D-Si_Ech"], mm["D-Si_Ych300"], fill_value=0, bounds_error=False)
        self.Y_D_Si_ch600 = interp1d(mm["D-Si_Ech"], mm["D-Si_Ych600"], fill_value=0, bounds_error=False)
        self.Y_D_C = interp1d(mm["D-C_E"], mm["D-C_Y"], fill_value=0, bounds_error=False)
        self.Y_D_C_ch = interp1d(mm["D-C_Ech"], mm["D-C_Ychsurf"], fill_value=0, bounds_error=False)
        self.Y_D_SiC_C = interp1d(mm["D-SiC,C_E"], mm["D-SiC,C_Y"], fill_value=0, bounds_error=False)
        self.Y_D_SiC_Cch = interp1d(mm["D-SiC,C_Ech"], mm["D-SiC,C_Ychsurf"], fill_value=0, bounds_error=False)
        self.Y_D_SiC_Si = interp1d(mm["D-SiC,Si_E"], mm["D-SiC,Si_Y"], fill_value=0, bounds_error=False)
        self.Y_D_SiC_Sich = interp1d(mm["D-SiC,Si_Ech"], mm["D-SiC,Si_Ychsurf"], fill_value=0, bounds_error=False)
        self.Y_C_Si = interp1d(mm["C-Si_E"], mm["C-Si_Y"], fill_value=0, bounds_error=False)
        self.Y_C_C = interp1d(mm["C-C_E"], mm["C-C_Y"], fill_value=0, bounds_error=False)
        self.Y_C_SiC_C = interp1d(mm["C-SiC,C_E"], mm["C-SiC,C_Y"], fill_value=0, bounds_error=False)
        self.Y_C_SiC_Si = interp1d(mm["C-SiC,Si_E"], mm["C-SiC,Si_Y"], fill_value=0, bounds_error=False)

        # Adding Si impact ion yields.
        self.Y_Si_C = interp1d(mm["Si-C_E"], mm["Si-C_Y"], fill_value=0, bounds_error=False)
        self.Y_Si_Si = interp1d(mm["Si-Si_E"], mm["Si-Si_Y"], fill_value=0, bounds_error=False)
        self.Y_Si_SiC_C = interp1d(mm["Si-SiC,C_E"], mm["Si-SiC,C_Y"], fill_value=0, bounds_error=False)
        self.Y_Si_SiC_Si = interp1d(mm["Si-SiC,Si_E"], mm["Si-SiC,Si_Y"], fill_value=0, bounds_error=False)

    def load_w_yields(self):
        """
        We add in the functionality to use similar C and D on W yields to calculate the sputtering flux of W.
        """

    def run_mm_model(self):

        # Find out what te is at the limiter location.
        te = self.fte(self.rlim)
        eimp = 3 * te + 2 * te * self.timult
        print("Te = {:.2f}   Eimp = {:.2f}".format(te, eimp))

        # Pick yields and assign to friendlier variable names.
        #A1 = self.Y_D_Si(eimp)
        A1 = self.Y_D_SiC_Si(eimp)
        A2 = self.Y_D_SiC_Sich(eimp)
        A3 = self.Y_C_SiC_Si(eimp)
        A4 = self.Y_D_SiC_C(eimp)
        A5 = self.Y_D_SiC_Cch(eimp)
        A6 = self.Y_C_SiC_C(eimp)
        A7 = self.Y_D_Si(eimp)
        A8 = self.Y_D_Si_ch25(eimp)
        A9 = self.Y_C_Si(eimp)
        A10 = self.Y_D_C(eimp)
        A11 = self.Y_D_C_ch(eimp)
        A12 = self.Y_C_C(eimp)

        # Adding Si impact yields.
        A13 = self.Y_Si_SiC_Si(eimp)
        A14 = self.Y_Si_SiC_C(eimp)
        A15 = self.Y_Si_Si(eimp)
        A16 = self.Y_Si_C(eimp)


        # Our system of 8 nonlinear equations.
        # Variables are assigned as such:
        # x1: Y_SiC, Si
        # x2: Y_SiC, C
        # x3: Y_Si
        # x4: Y_C
        # x5: Y_Si, tot
        # x6: Y_C, tot
        # x7: conc_C
        # x8: conc_Si
        # x9: fC (=x6)
        # x10: fSi (=x5)
        R = 0.1  # Apparently fSi = fC * R
        def equations(vars):
            x1, x2, x3, x4, x5, x6, x7, x8 = vars
            #if include_si_sput:
            eqs = [
                A1 + A2 + x6 * A3 + x5 * A13 - x1,
                A4 + A5 + x6 * A6 + x5 * A14 - x2,
                A7 + A8 + x6 * A9 + x5 * A15 - x3,     # Fixed a typo here, had * instead of +.
                A10 + A11 + x6 * A12 + x5 * A16 - x4,  # Fixed a typo here, had * instead of +.
                (1 - x7 - x8) * x1 + x8 * x3 - x5,
                (1 - x7 - x8) * x2 + x7 * x4 - x6,
                (1 - R) * x6 / x4 - x7,
                (1 - (1 - R) * x6 / x4) * (x2 - x1) / (x3 + x2 - x1) - x8]
            #else:
            #    eqs = [
            #        A1 + A2 + x6 * A3 - x1,
            #        A4 + A5 + x6 * A6 - x2,
            #        A7 + A8 * x6 * A9 - x3,
            #        A10 + A11 * x6 * A12 - x4,
            #        (1 - x7 - x8) * x1 + x8 * x3 - x5,
            #        (1 - x7 - x8) * x2 + x7 * x4 - x6,
            #        (1 - R) * x6 / x4 - x7,
            #        (1 - (1 - R) * x6 / x4) * (x2 - x1) / (x3 + x2 - x1) - x8]
            return eqs

        # Solve the system of equations and pull out fC.
        #guess = [0.05, 0.05, 0.05, 0.05, 0.002, 0.05, 0.5, 0.5]
        guess = [0.0005, 0.015, 0.02, 0.02, 0.002, 0.02, 0.7, 0.08]
        xs = fsolve(equations, guess)

        for i in range(0, len(xs)):
            print("x{}: {}".format(i+1, xs[i]))

        fc_sic = xs[5]
        fsi_sic = xs[4]

        # For graphite comparisons it's a simple equation.
        fc_gph = (A10 + A11) / (1 - A12)

        # Calculate a rough Zeff.
        zeff_sic = xs[5] * 6 + xs[4] * 14 + (1 - xs[5] - xs[4]) * 1
        zeff_gph = fc_gph * 6 + (1-fc_gph) * 1

        # Store in class.
        self.fc_sic = fc_sic
        self.fsi_sic = fsi_sic
        self.fc_gph = fc_gph
        self.zeff_sic = zeff_sic
        self.zeff_gph = zeff_gph
        self.Y_Si_tot = xs[6]
        self.Y_C_tot = xs[7]
        self.Y_C = xs[3]

        # Return results as a dictionary.
        return {"fc_sic":fc_sic, "fsi_sic":fsi_sic, "fc_gph":fc_gph,
            "zeff_sic":zeff_sic, "zeff_gph":zeff_gph, "Y_Si_tot":xs[6],
            "Y_C_tot":xs[7], "Y_C":xs[3]}

    def calc_neut_flux(self, mat):

        # Grab corresponding yield.
        if mat == "c_sic":
            fz = self.fc_sic
        elif mat == "si_sic":
            fz = self.fsi_sic
        elif mat == "c":
            fz = self.fc_gph

    def run_engelhardt(self, mat, dperp=[1, 1, 5, 10],
        conns=[100, 50, 10, 5], vels=[500, 500, 500, 500],
        match_peter=False, mode="diffusive"):
        """

        vels (list, float): Inputs for the radial velocities in each region.
        match_peter (bool): One time thing to match Peter's Excel file.
        mode (str): Either "diffusive" or "convective".
        """

        # Pull out some items for easy access.
        r0 = self.r0
        rsep = 0.0
        rlim = self.rlim
        riz = self.riz
        rwall = self.rwall
        ne = self.ne
        te = self.te
        ti = self.ti

        # Choose from either C (SiC), Si (SiC), or just graphite. Additionally,
        # key insight is that flux_neut = fz * ne!
        if mat == "c_sic":
            fz = self.fc_sic
        elif mat == "si_sic":
            fz = self.fsi_sic
        elif mat == "c":
            fz = self.fc_gph
        neut_flux = fz * ne[3]
        print("neut_flux = {:.2e}".format(neut_flux))

        # Sound speed.
        mi = 2 * 931.49e6
        cs = np.sqrt((te+ti)/mi) * 3e8

        # Root for the equations.
        conns = np.array(conns)
        dperp = np.array(dperp)
        x = np.sqrt(cs / (dperp * conns))

        if match_peter:
            x = [1, 1/0.03, 1/0.03, 1/0.01]

        nr = 50
        rsd = np.linspace(rlim, rwall, nr)
        rsc = np.linspace(riz, rlim, nr)
        rsb = np.linspace(rsep, riz, nr)
        rsa = np.linspace(r0, rsep)

        if mode == "diffusive":

            # Region d solution.
            # Hard boundary condition.
            #nzd = fz * ne[3] * (np.exp(x[3]*rsd) - np.exp(2*x[3]*rwall) * np.exp(-x[3]*rsd)) \
            #    / (np.exp(x[3]*rlim) - np.exp(2*x[3]*rwall) * np.exp(-x[3]*rlim))

            # Soft boundary condition (assume normal decaying exponential).
            nzd = fz * ne[3] * np.exp(-x[3]*(rsd-rlim))

            # Region c solution.
            C1 = (fz * ne[3] * np.exp(x[2]*rlim) * np.exp(-x[2]*riz) - neut_flux / dperp[2]) \
                / (x[2]*np.exp(x[2]*riz) + x[2]*np.exp(2*x[2]*rlim) * np.exp(-x[2]*riz))
            nzc = C1 * (np.exp(x[2]*rsc) - np.exp(2*x[2]*rlim) * np.exp(-x[2]*rsc)) \
                + fz * ne[3] * np.exp(x[2]*rlim) * np.exp(-x[2]*rsc)

            # Region b solution.
            C2 = (C1 * (np.exp(x[2]*riz) - np.exp(2*x[2]*rlim) * np.exp(-x[2]*riz)) \
                + fz * ne[3] * np.exp(x[2]*rlim) * np.exp(-x[2]*riz)) \
                / (np.exp(x[1]*riz) + np.exp(2*x[2]*rsep) * np.exp(-x[1]*riz))
            nzb = C2 * (np.exp(x[1]*rsb) + np.exp(2*x[1]*rsep) * np.exp(-x[1]*rsb))

        elif mode == "convective":

            lamb = conns * vels / cs

            # Region d solution.
            nzd = fz * ne[3] * np.exp((rlim-rsd) / lamb[3])

            # Region c solution.
            nzc = fz * ne[3] * np.exp((rlim-rsc) / lamb[2])

            # Region b solution.
            nzb = fz * ne[3] * np.exp((rlim-riz) / lamb[2] + (riz - rsb) / lamb[1])


        # Region a solution.

        nza = np.full(nr, nzb[0])

        print("Core {} density: {:.2e} m-3".format(mat, nza[0]))

        # Zip into one full profile.
        rs = np.append(rsa, np.append(rsb, np.append(rsc, rsd)))
        nz = np.append(nza, np.append(nzb, np.append(nzc, nzd)))
        self.rs = rs
        self.nz = nz

        return {"rs":rs, "nz":nz}

    def calc_zeff_prof_sic(self, rs, nz_c, nz_si):

        # First get the corresponding ne value at each nz location, then the
        # fractions of each ion.
        ne = self.fne(rs)
        fc = nz_c / ne
        fsi = nz_si / ne

        zeff_sic = fc * 6 + fsi * 14 + (1 - fsi - fc) * 1
        return zeff_sic

    def calc_zeff_prof_gph(self, rs, nz_c):
        ne = self.fne(rs)
        fc = nz_c / ne
        zeff_gph = fc * 6 + (1-fc) * 1
        return zeff_gph
