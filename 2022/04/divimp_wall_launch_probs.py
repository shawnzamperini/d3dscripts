# This script loads in some data from the 167196 DIVIMP grid to find the launch
# probabilities for C and Si from graphite and SiC walls. The output is passed
# into input option N16.
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import oedge_plots


# Option to load fC impinging upon the wall elements from a previous run. The
# idea is to iterate this multiple times and hopefully converge on a solution.
# Generally set to True, as False just uses a W case to load the geometry of the
# situation and starts with a constant fC = 0.02 at every surface.
load_fc = False
max_fC = 0.05  # Only applies if load_fc = True
fill_fC = 0.02
impact_charge = 2.0
shelf_only = True  # Only print out data for the shelf tiles (indices 196-256)

if shelf_only:
    print("Only printing out data for the shelf tiles")

if load_fc:
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-grap.nc"
else:
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mrc-shifted-nodrift-2.nc"

# Can extract all the needed info from the netCDF file. See notes for details on
# which row contains what data.
op = oedge_plots.OedgePlots(ncpath)
wallpt = op.nc.variables["WALLPT"][:]

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

refl = 0.1
def calc_yield(te, atom, ti_mult=2.0, fC=0.02, cZ=1):
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


# Pull out the data we actually want from the wall.
tes = wallpt[28]
keep = np.nonzero(tes)
tes = tes[keep]
tis = wallpt[29][keep]
nes = wallpt[30][keep]
lens = wallpt[6][keep]
rs = wallpt[0][keep]

# Flux just the density times cs.
cs = np.sqrt((tes + tis) * 1.609e-19 / (2.0 * 1.66e-27))
fluxes = nes * cs

if load_fc:

    # Unfortunately fC isn't easily accessed from WALLPT, as the data is only
    # available on target segments, i.e. where the grid exists. What we can do
    # it get the data for wall elements that correspond to a target element and
    # assign fC, and then for all other segments just assign fC to the lowest
    # nonzero value found elsewhere. Unfortunately the IK, IR indices supposedly
    # coded into wallpt[25] and wallpt[26] are incomplete, so we take a round
    # about way to get the flux at each target element.
    wall_targ_idxs = wallpt[17][keep]

    # Ring, knot indices at each target element.
    irds = op.nc.variables["IRDS"][:]
    irks = op.nc.variables["IKDS"][:]
    ddlims = op.nc.variables["DDLIMS"][:]
    nrs = int(op.nc.variables["NRS"][:])
    nks = op.nc.variables["NKS"][:]
    cion = int(op.nc.variables["CION"][:])

    # Average charge of the impurities at a location is the average weighted by
    # each state's respective density.
    #print("Calculating average charges...")
    #avg_charges = np.ones(ddlims[0].shape)
    #for ir in range(0, nrs):
    #    for ik in range(0, int(nks[ir])):
    #        avg_charge = 0
    #        tot_dens = ddlims[:,ir,ik].sum()
    #        if tot_dens == 0:
    #            continue
    #        for charge in range(1, cion+1):
    #
    #            # Reminder index 0 is neutrals.
    #            avg_charge += charge * ddlims[charge][ir][ik] / tot_dens
    #        avg_charges[ir][ik] = avg_charge

    avg_charges = np.zeros(len(wall_targ_idxs))
    tot_dens = np.zeros(len(wall_targ_idxs))
    fC = np.zeros(len(wall_targ_idxs))

    # i here is still the wall_idx (zero indexed though, so +1 for the real index).
    for i in range(0, len(wall_targ_idxs)):
        targ = int(wall_targ_idxs[i]) - 1  # Go from 1 to 0 indexed.
        ring = irds[targ]
        knot = irks[targ]

        # With the ring, knot index at the wall/target location, we now can
        # easily access anything. First we can grab the average charge state
        # (as a float!) impinging at this location. For DDLIMS I think the first
        # two "charges" are primary and then total neutrals, and then last entries
        # are for each charge state.
        avg_charge = 0.0
        tot_dens[i] = ddlims[2:,ring,knot].sum()
        if tot_dens[i] == 0:
            continue
        for charge in range(1, cion+1):
            avg_charge += charge * ddlims[charge][ring][knot] / tot_dens[i]
        avg_charges[i] = avg_charge

    # Where the average charge = 0, just assign the minimum nonzero average
    # charge calculated.
    avg_charges[avg_charges==0] = avg_charges[avg_charges!=0].min()

    # Fraction of carbon simply the total C density (apply absfac) divided by ne
    fC = tot_dens * op.absfac / nes
    fC = np.clip(fC, 0.0, max_fC)

    # For area where fC = 0.0, i.e. no corresponding target to match to, assign
    # just the mean exlucindg zeros.
    fC[fC==0] = fC[fC==0].mean()

else:
    fC = np.full(len(tes), fill_fC)
    avg_charges = np.full(len(tes), impact_charge)

ti_mult = 1.0
grap_weighted = np.zeros(len(tes))
sic_c_weighted = np.zeros(len(tes))
sic_si_weighted = np.zeros(len(tes))
grap_ys = np.zeros(len(tes))
sic_c_ys = np.zeros(len(tes))
sic_si_ys = np.zeros(len(tes))

# Calculate the mixed material yields.
for i in range(0, len(tes)):

    # i is the wall_idx, so if we only want the shelf points we can do that here.
    if shelf_only:
        if i < 196 or i > 256:
            continue

    grap_ys[i] = calc_yield(tes[i], "Conly", ti_mult=ti_mult, fC=fC[i], cZ=avg_charges[i])
    sic_c_ys[i] = calc_yield(tes[i], "C", ti_mult=ti_mult, fC=fC[i], cZ=avg_charges[i])
    sic_si_ys[i] = calc_yield(tes[i], "Si", ti_mult=ti_mult, fC=fC[i], cZ=avg_charges[i])
    grap_flux = grap_ys[i] * fluxes[i]
    sic_c_flux = sic_c_ys[i] * fluxes[i]
    sic_si_flux = sic_si_ys[i] * fluxes[i]

    # Weigh by the lengths. Note here we are not actually changing units,
    # more properly we normalize the lengths and multiply but since we normalize
    # later it doesn't matter. Further multiply by 2piR so that we correctly
    # weight by the full area of each segment. Units at the end are conceptually
    # of atoms/tor-m/s.
    grap_weighted[i] = grap_flux * lens[i] * 2 * np.pi * rs[i]
    sic_c_weighted[i] = sic_c_flux * lens[i] * 2 * np.pi * rs[i]
    sic_si_weighted[i] = sic_si_flux * lens[i] * 2 * np.pi * rs[i]

# ABSFAC will be the maximum value of the weighted values. Units of atoms/tor-m/s.
absfac_grap = grap_weighted.max()
absfac_sic_c = sic_c_weighted.max()
absfac_sic_si = sic_si_weighted.max()

# To do probabilities just divide by the sum.
grap_prob = grap_weighted / grap_weighted.sum()
sic_c_prob = sic_c_weighted / sic_c_weighted.sum()
sic_si_prob = sic_si_weighted / sic_si_weighted.sum()

# Save to an Excel file.
wall_idx = np.arange(1, len(tes)+1)
df = pd.DataFrame({"wall_idx1":wall_idx, "wall_idx2":wall_idx, "tes":tes, "nes":nes,
    "grap_prob":grap_prob, "sic_c_prob":sic_c_prob, "sic_si_prob":sic_si_prob,
    "absfac_grap":absfac_grap, "absfac_sic_c":absfac_sic_c,
    "absfac_sic_si":absfac_sic_si})
df.to_excel("/Users/zamperini/My Drive/Research/Documents/2022/04/wall_launch_probs.xlsx")
