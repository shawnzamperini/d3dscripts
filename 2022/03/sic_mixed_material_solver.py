# This script solves a system of nonlinear equations for Tyler's SiC mixed
# material model under the assumption that the incoming carbon flux is equal
# to the outgoing carbon flux. This mathematically means fC = Y_C_tot.
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt


# Constants to determine the constituent yields that we already know.
tes = np.linspace(1, 60, 100)
ti_mult = 1.0
R = 0.1
guess = [0.05, 0.05, 0.05, 0.05, 0.002, 0.05, 0.5, 0.5]
graphite_only = False
include_si_sput = False

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
Y_D_Si = interp1d(mm["D-Si_E"], mm["D-Si_Y"], fill_value=0, bounds_error=False)
Y_D_Si_ch25 = interp1d(mm["D-Si_Ech"], mm["D-Si_Ych25"], fill_value=0, bounds_error=False)
Y_D_Si_ch300 = interp1d(mm["D-Si_Ech"], mm["D-Si_Ych300"], fill_value=0, bounds_error=False)
Y_D_Si_ch600 = interp1d(mm["D-Si_Ech"], mm["D-Si_Ych600"], fill_value=0, bounds_error=False)
Y_D_C = interp1d(mm["D-C_E"], mm["D-C_Y"], fill_value=0, bounds_error=False)
Y_D_C_ch = interp1d(mm["D-C_Ech"], mm["D-C_Ychsurf"], fill_value=0, bounds_error=False)
Y_D_SiC_C = interp1d(mm["D-SiC,C_E"], mm["D-SiC,C_Y"], fill_value=0, bounds_error=False)
Y_D_SiC_Cch = interp1d(mm["D-SiC,C_Ech"], mm["D-SiC,C_Ychsurf"], fill_value=0, bounds_error=False)
Y_D_SiC_Si = interp1d(mm["D-SiC,Si_E"], mm["D-SiC,Si_Y"], fill_value=0, bounds_error=False)
Y_D_SiC_Sich = interp1d(mm["D-SiC,Si_Ech"], mm["D-SiC,Si_Ychsurf"], fill_value=0, bounds_error=False)
Y_C_Si = interp1d(mm["C-Si_E"], mm["C-Si_Y"], fill_value=0, bounds_error=False)
Y_C_C = interp1d(mm["C-C_E"], mm["C-C_Y"], fill_value=0, bounds_error=False)
Y_C_SiC_C = interp1d(mm["C-SiC,C_E"], mm["C-SiC,C_Y"], fill_value=0, bounds_error=False)
Y_C_SiC_Si = interp1d(mm["C-SiC,Si_E"], mm["C-SiC,Si_Y"], fill_value=0, bounds_error=False)

# Adding Si impact ion yields.
Y_Si_C = interp1d(mm["Si-C_E"], mm["Si-C_Y"], fill_value=0, bounds_error=False)
Y_Si_Si = interp1d(mm["Si-Si_E"], mm["Si-Si_Y"], fill_value=0, bounds_error=False)
Y_Si_SiC_C = interp1d(mm["Si-SiC,C_E"], mm["Si-SiC,C_Y"], fill_value=0, bounds_error=False)
Y_Si_SiC_Si = interp1d(mm["Si-SiC,Si_E"], mm["Si-SiC,Si_Y"], fill_value=0, bounds_error=False)


# Go through one Te at a time.
fcs_sic = np.zeros(len(tes))
fsis_sic = np.zeros(len(tes))
fcs_gph = np.zeros(len(tes))
zeff_sic = np.zeros(len(tes))
zeff_gph = np.zeros(len(tes))
for i in range(0, len(tes)):
    print("Te = {:.2f}".format(tes[i]))
    eimp = 3 * tes[i] + 2 * tes[i] * ti_mult

    # Pick yields and assign to friendlier variable names.
    A1 = Y_D_Si(eimp)
    A2 = Y_D_SiC_Sich(eimp)
    A3 = Y_C_SiC_Si(eimp)
    A4 = Y_D_SiC_C(eimp)
    A5 = Y_D_SiC_Cch(eimp)
    A6 = Y_C_SiC_C(eimp)
    A7 = Y_D_Si(eimp)
    A8 = Y_D_Si_ch25(eimp)
    A9 = Y_C_Si(eimp)
    A10 = Y_D_C(eimp)
    A11 = Y_D_C_ch(eimp)
    A12 = Y_C_C(eimp)

    # Adding Si impact yields.
    A13 = Y_Si_SiC_Si(eimp)
    A14 = Y_Si_SiC_C(eimp)
    A15 = Y_Si_Si(eimp)
    A16 = Y_Si_C(eimp)


    # Our system of 8 nonlinear equations.
    # Variables are assigned as such:
    # ...
    def equations(vars):
        x1, x2, x3, x4, x5, x6, x7, x8 = vars
        if include_si_sput:
            eqs = [
                A1 + A2 + x6 * A3 + x5 * A13 - x1,
                A4 + A5 + x6 * A6 + x5 * A14 - x2,
                A7 + A8 * x6 * A9 + x5 * A15 - x3,
                A10 + A11 * x6 * A12 + x5 * A16 - x4,
                (1 - x7 - x8) * x1 + x8 * x3 - x5,
                (1 - x7 - x8) * x2 + x7 * x4 - x6,
                (1 - R) * x6 / x4 - x7,
                (1 - (1 - R) * x6 / x4) * (x2 - x1) / (x3 + x2 - x1) - x8]
        else:
            eqs = [
                A1 + A2 + x6 * A3 - x1,
                A4 + A5 + x6 * A6 - x2,
                A7 + A8 * x6 * A9 - x3,
                A10 + A11 * x6 * A12 - x4,
                (1 - x7 - x8) * x1 + x8 * x3 - x5,
                (1 - x7 - x8) * x2 + x7 * x4 - x6,
                (1 - R) * x6 / x4 - x7,
                (1 - (1 - R) * x6 / x4) * (x2 - x1) / (x3 + x2 - x1) - x8]
        return eqs

    # Solve the system of equations and pull out fC.
    xs = fsolve(equations, guess)
    fcs_sic[i] = xs[5]
    fsis_sic[i] = xs[4]

    # For graphite comparisons it's a simple equation.
    fcs_gph[i] = (A10 + A11) / (1 - A12)

    # Calculate a rough Zeff.
    zeff_sic[i] = xs[5] * 6 + xs[4] * 14 + (1 - xs[5] - xs[4]) * 1
    zeff_gph[i] = fcs_gph[i] * 6 + (1-fcs_gph[i]) * 1

# Remove anything negative (solver failed).
mask = fcs_sic > 0
tes = tes[mask]
fcs_sic = fcs_sic[mask]
fsis_sic = fsis_sic[mask]
fcs_gph = fcs_gph[mask]
zeff_sic = zeff_sic[mask]
zeff_gph = zeff_gph[mask]

# For SiC 7eV seems to be the cutoff.
sic_mask = tes > 6.8

plt.rcParams['font.family'] = 'Century Gothic'
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13,4))
ax1.plot(tes[sic_mask], fcs_sic[sic_mask]*100, label="SiC", lw=3, color="tab:purple")
ax1.plot(tes, fcs_gph*100, label="Graphite", lw=3, color="tab:red")
ax1.grid()
ax1.legend(fontsize=16)
ax1.tick_params(labelsize=14)
ax1.set_ylim([0, 10])
ax1.set_xlabel("Te (eV)", fontsize=16)
ax1.set_ylabel("Fraction of C in plasma (%)", fontsize=16)
ax2.plot(tes[sic_mask], fsis_sic[sic_mask]*100, label="SiC", lw=3, color="tab:purple")
ax2.grid()
ax2.legend(fontsize=16)
ax2.tick_params(labelsize=14)
ax2.set_ylim([0, 10])
ax2.set_xlabel("Te (eV)", fontsize=16)
ax2.set_ylabel("Fraction of Si in plasma (%)", fontsize=16)
ax3.plot(tes[sic_mask], zeff_gph[sic_mask], label="Graphite", lw=3, color="tab:red")
ax3.plot(tes[sic_mask], zeff_sic[sic_mask], label="SiC", lw=3, color="tab:purple")
ax3.grid()
ax3.legend(fontsize=16)
ax3.tick_params(labelsize=14)
ax3.set_ylim([1, 2])
ax3.set_xlabel("Te (eV)", fontsize=16)
ax3.set_ylabel("Zeff", fontsize=16)
fig.tight_layout()
fig.show()
