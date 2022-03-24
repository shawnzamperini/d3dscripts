# This script does a simple implemetation of the SiC mixed material model.
# The goal is a plot of sputtered particles vs. Te at the limiter surface.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D


# Inputs.
refl = 0.1
cZ = 2.0

# Load in mixed material data.
mmpath = "/Users/zamperini/My Drive/Research/Documents/2022/02/mixed_material.xlsx"
columns = ["D-Si_E", "D-Si_Y", "D-Si_Ech", "D-Si_Ych25", "D-Si_Ych300",
  "D-Si_Ych600", "D-C_E", "D-C_Y", "D-C_Ech", "D-C_Ychsurf", "D-C_Ychphys",
  "D-SiC,C_E", "D-SiC,C_Y", "D-SiC,C_Ech", "D-SiC,C_Ychsurf", "D-SiC,Si_E",
  "D-SiC,Si_Y", "D-SiC,Si_Ech", "D-SiC,Si_Ychsurf", "C-Si_E",
  "C-Si_Y", "C-C_E", "C-C_Y", "C-SiC,C_E", "C-SiC,C_Y", "C-SiC,Si_E",
  "C-SiC,Si_Y"]
mm = pd.read_excel(mmpath, skiprows=6, header=None, names=columns,
    usecols="A:F,H:L,N:Q,S:V,X,Y,AA,AB,AD,AE,AG,AH")

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

def calc_yield(te, atom, ti_mult=2.0, fC=0.02):
    """
    Calculate yield using mixed material model. atom is one of "C" or "Si" for
    SiC, or simple "Conly" or "Sionly" if you just want the monoatomic yields.
    """

    # Impact energies of deuterium and carbon.
    ED = 3 * te + 2 * ti_mult * te
    EC = 3 * cZ * te + 2 * ti_mult * te

    # Yield calculations from Abrams NF 2021.
    Y_C = Y_D_C(ED) + Y_D_C_ch(ED) + fC * Y_C_C(EC)
    Y_Si = Y_D_Si(ED) + Y_D_Si_ch25(ED) + fC * Y_C_Si(EC)
    Y_SiC_C = Y_D_SiC_C(ED) + Y_D_SiC_Cch(ED) + fC * Y_C_SiC_C(EC)
    Y_SiC_Si = Y_D_SiC_Si(ED) + Y_D_SiC_Sich(ED) + fC * Y_C_SiC_Si(EC)

    # Surface concentrations from Abrams NF 2021.
    if Y_C == 0.0:
        conc_C = 0.0
    else:
        conc_C = (1 - refl) * fC / Y_C
    if (Y_Si + Y_SiC_C - Y_SiC_Si) == 0.0:
        conc_Si = 0.0
    else:
        conc_Si = (1 - conc_C) * (Y_SiC_C - Y_SiC_Si) / (Y_Si + Y_SiC_C - Y_SiC_Si)
        if conc_Si < 0:
            conc_Si = 0
    conc_SiC = 1 - conc_C - conc_Si

    #print("{:.1f} {:.1f} {:.1f}: {:.3f} {:.3f} {:.3f}".format(te, ED, EC, conc_C, conc_Si, Y_C))
    #print("{:.1f} {:.1f} {:.1f}: {:.3f} {:.3f} {:.3f}".format(te, ED, EC, conc_C, Y_C, Y_D_C(ED)))

    # Total yields.
    Y_Ctot = conc_SiC * Y_SiC_C + conc_C * Y_C
    Y_Sitot = conc_SiC * Y_SiC_Si + conc_Si * Y_Si

    if atom.lower() == "c":
        return Y_Ctot
    elif atom.lower() == "si":
        return Y_Sitot
    elif atom.lower() == "conly":
        return Y_C
    elif atom.lower() == "sionly":
        return Y_Si

# Yields for SiC at various Ti multipliers.
tes = np.linspace(0, 50, 100)
nscan = 10
mults = np.linspace(1, 5, nscan)
Ys_SiC_C = np.zeros((len(mults), len(tes)))
Ys_SiC_Si = np.zeros((len(mults), len(tes)))
Ys_C = np.zeros((len(mults), len(tes)))
Ys_Si = np.zeros((len(mults), len(tes)))
for j in range(0, len(mults)):
    for i in range(0, len(tes)):
        Ys_SiC_C[j][i] = calc_yield(tes[i], "C", ti_mult=mults[j])
        Ys_SiC_Si[j][i] = calc_yield(tes[i], "Si", ti_mult=mults[j])
        Ys_C[j][i] = calc_yield(tes[i], "Conly", ti_mult=mults[j])
        Ys_Si[j][i] = calc_yield(tes[i], "Sionly", ti_mult=mults[j])
Ys_SiC_C = np.ma.masked_where(Ys_SiC_C<=0, Ys_SiC_C)
Ys_SiC_Si = np.ma.masked_where(Ys_SiC_Si<=0, Ys_SiC_Si)
Ys_C = np.ma.masked_where(Ys_C<=0, Ys_C)
Ys_Si = np.ma.masked_where(Ys_Si<=0, Ys_Si)

# Yields for different fractions of carbon.
cfracs = np.linspace(0, 0.05, 10)
Ys_SiC_C2 = np.zeros((len(cfracs), len(tes)))
Ys_SiC_Si2 = np.zeros((len(cfracs), len(tes)))
Ys_C2 = np.zeros((len(cfracs), len(tes)))
Ys_Si2 = np.zeros((len(cfracs), len(tes)))
for j in range(0, len(cfracs)):
    for i in range(0, len(tes)):
        Ys_SiC_C2[j][i] = calc_yield(tes[i], "C", fC=cfracs[j])
        Ys_SiC_Si2[j][i] = calc_yield(tes[i], "Si", fC=cfracs[j])
        Ys_C2[j][i] = calc_yield(tes[i], "Conly", fC=cfracs[j])
        Ys_Si2[j][i] = calc_yield(tes[i], "Sionly", fC=cfracs[j])
Ys_SiC_C2 = np.ma.masked_where(Ys_SiC_C2<=0, Ys_SiC_C2)
Ys_SiC_Si2 = np.ma.masked_where(Ys_SiC_Si2<=0, Ys_SiC_Si2)
Ys_C2 = np.ma.masked_where(Ys_C2<=0, Ys_C2)
Ys_Si2 = np.ma.masked_where(Ys_Si2<=0, Ys_Si2)

# Colormap.
plt.rcParams['font.family'] = 'Century Gothic'
cmap = plt.get_cmap('inferno')
colors = cmap(np.linspace(0, 0.9, len(mults)))
custom_lines = [Line2D([0], [0], color=colors[0], linestyle="--"),
    Line2D([0], [0], color=colors[0])]
fontsize = 14

# Plot of yields vs Te for different conditions.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
for j in range(0, nscan):
    ax1.plot(tes, Ys_SiC_C[j], color=colors[j])
    ax1.plot(tes, Ys_C[j], color=colors[j], linestyle="--")
    ax2.plot(tes, Ys_SiC_C2[j], color=colors[j])
    ax2.plot(tes, Ys_C2[j], color=colors[j], linestyle="--")
#ax1.set_yscale("log")
ax1.set_xlabel("Te (eV)", fontsize=fontsize)
ax1.set_ylabel("Yield", fontsize=fontsize)
ax1.grid(axis="x", which="major")
ax1.grid(axis="y", which="both")
ax1.set_title(r"$\mathdefault{C^{2+}\ f_C}$ = 2%", fontsize=fontsize)
ax2.set_xlabel("Te (eV)", fontsize=fontsize)
ax2.grid(axis="x", which="major")
ax2.grid(axis="y", which="both")
ax2.set_title(r"$\mathdefault{C^{2+}\ T_i\ =\ 2T_e}$", fontsize=fontsize)
ax1.set_ylim(0.0, 0.08)
ax2.set_ylim(0.0, 0.08)
ax2.set_yticklabels([])
ax1.legend(custom_lines, ["Graphite", "SiC"], loc="upper left", fontsize=fontsize)
fig.tight_layout()
fig.show()
