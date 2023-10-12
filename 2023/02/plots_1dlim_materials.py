import LimPlots
import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append("../../2022/08")
import EngelhardtModel

mi = 2 * 931.49e6

def trim(lp, Z):
    xouts = lp.nc['XOUTS'][:].data
    youts = lp.nc['YOUTS'][:].data
    xwids = lp.nc["XWIDS"][:].data
    mask = xwids != 0
    xouts = xouts[mask]
    xwids = xwids[mask]
    xkeep_min = np.nonzero(xouts)[0].min()
    xkeep_max = np.nonzero(xouts)[0].max()
    ykeep_min = np.nonzero(youts)[0].min()
    ykeep_max = np.nonzero(youts)[0].max()
    xouts = xouts[xkeep_min:xkeep_max]
    youts = youts[ykeep_min:ykeep_max]
    xwids = xwids[xkeep_min:xkeep_max]
    Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
    yabsorb1a = float(lp.nc["yabsorb1a"][:].data)
    yabsorb2a = float(lp.nc["yabsorb2a"][:].data)
    ykeep = np.where(np.logical_and(youts >= yabsorb2a, youts <= yabsorb1a))[0]
    youts = youts[ykeep]
    Z = Z[ykeep, :]
    return Z


def get_zeff(lp, ion, lambda_ne, lambda_te, nesep=1e19, tesep=100, sic=True, ftot2_1=0, ftot2_2=0, fytot2_1=0,
             fytot2_2=0, riz=0.01, rlim=0.10, rwall=0.12):

    try:
        nz_dict = lp.plot_par_rad("nz", 20, charge="all", showplot=False)
        nz = nz_dict["Z"].data[0]
    except IndexError:
        return {"rmrs": [0], "zeff": [0], "zeff_gph": [0]}


    # Convert to R-Rsep.
    r = 0.10 - nz_dict["Y"][0]

    # We use just the mixed-material portion of the EngelhardtModel that I coded up to grab what the fraction of C and
    # Si at the surface is, which we then use to scale the density results to agree with the values at the surface.
    # This is because the value of 1DLIM a) doesn't self-consistently solve for fC and fSi and b) the ABSFAC seems wrong
    # anyways.
    # if sic:
    em = EngelhardtModel.EngelhardtModel()
    em.set_geometry(riz=riz, rlim=rlim, rwall=rwall)
    em.set_ne(nesep=nesep, lambda_ne=lambda_ne)
    em.set_te_ti(tesep=tesep, lambda_te=lambda_te, timult=3.0)
    em.load_mm_model()
    out1 = em.run_mm_model()

    # Determine the appropriate scaling factor to bring into agreement with the mixed-material model.
    lim_idx = np.where(np.abs(r - rlim) == np.abs(r - rlim).min())
    lp_fact = nz[lim_idx[0][0]]

    if ion == "C":
        maxz = 7
        if sic:
            fact = out1["fc_sic"] * em.fne(rlim)
    elif ion == "Si":
        maxz = 15
        if sic:
            fact = out1["fsi_sic"] * em.fne(rlim)
    elif ion == "W":
        maxz = 74
        te = em.fte(rlim)
        fact = calc_yeff(te, te, 1) * em.fne(rlim)
        absfac = 1.0

    # if not sic:
    #     fact = 1.0
    #     lp_fact = 1.0

    # lp_nz = nz * fact / lp_fact

    if ion in ["C", "Si"]:
        absfac = float(lp.nc["ABSFAC"][:].data)

    # W is a tricky thing since we are simulating it assuming a 2% or whatever amount background C. Only the C causes
    # any sputtering (secondary yields in the code with sputter option 2). But ABSFAC is only calculated from the
    # deuterium (primary) yields. The data is in the dat file though to calculate an ABSFAC from the carbon values,
    # and so we do that here.
    # GTOT2 = 0.5 * (FTOT2(1) + FTOT2(2))  This is, the average of the integrated flux values for each side.
    # GYTOT2 = 0.5 * (FYTOT2(1) + FYTOT2(2))   That is, the average of the flux*yield values for each side.
    # YEFF = GYTOT2 / GTOT2
    # ABSFAC = 2.0 * GTOT2 * YEFF
    # elif ion == "W":
    #     gtot2 = 0.5 * (ftot2_1 + ftot2_2)
    #     gytot2 = 0.5 * (fytot2_1 + fytot2_2)
    #     yeff = gytot2 / gtot2
    #     absfac = 2.0 * gtot2 * yeff

    ddlim3 = lp.nc["DDLIM3"][20, :, :, :].data * absfac * fact / lp_fact
    # ddlim3_all = trim(lp, ddlim3[1:].sum(axis=0)).data
    ne = lp.nc["CRNBS"][:].data
    ne = trim(lp, ne)[0]

    for charge in range(1, maxz):

        # Get the (unscaled) density for just this charge, trim it to just the
        # simulation volume.
        nz_charge = ddlim3[charge + 1]
        nz_charge = trim(lp, nz_charge)[0]

        # nz_dict = lp.plot_par_rad("nz", 20, charge=charge, showplot=False)
        # nz_charge = nz_dict["Z"].data[0]

        # Represent the density as what fraction of the total density it is, and
        # then multiply by our scaled density lp_nc.
        # charge_frac = nz_charge / ddlim3_all
        # nz_charge = nz_charge * lp_fact

        if charge == 1:
            # zeff = np.ones(nz_charge.shape)
            zeff = np.zeros(nz_charge.shape)
        zeff += charge ** 2 * nz_charge / ne

    # If we have the carbon run up, we can estimate the graphite contribution with just the graphite fraction.
    if sic and ion == "C":
        fact2 = out1["fc_gph"] * em.fne(rlim)
        return {"rmrs": r, "zeff": zeff, "zeff_gph":zeff * fact2 / fact}
    else:
        return {"rmrs": r, "zeff": zeff}

def calc_yeff(te, ti, zc):
    """
    It seems like the 1996 Eckstein data may have C-->W yield data, so we use the fit values from that to calculate
    the yield at the surface (yeff) and return that. This is used in calculating ABSFAC.
    """
    eth = 41.20
    etf = 66517.0
    q = 1.02
    ebd = 8.68
    energy = 2 * ti + 3 * te * zc

    x1 = energy / etf
    x2 = eth / energy
    yeff = q * (0.5 * np.log(1.0 + 1.2288 * x1)) / (x1 + 0.1728 * np.sqrt(x1) + 0.008 * x1**0.1504) * (1 - x2) \
           * (1 - x2) * (1.0 - x2**(2 / 3))

    return max(0, yeff)

# SiC (C) lambda = 5 cm
lp_sic_c_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb5.nc")
sic_c_5 = get_zeff(lp_sic_c_5, "C", 0.05, 0.05)

# SiC (Si) lambda = 5 cm
lp_sic_si_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb5.nc")
sic_si_5 = get_zeff(lp_sic_c_5, "Si", 0.05, 0.05)

# And likewise.
lp_sic_c_4 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb4.nc")
sic_c_4 = get_zeff(lp_sic_c_4, "C", 0.04, 0.04)
lp_sic_si_4 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb4.nc")
sic_si_4 = get_zeff(lp_sic_si_4, "Si", 0.04, 0.04)
lp_sic_c_3 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb3.nc")
sic_c_3 = get_zeff(lp_sic_c_3, "C", 0.03, 0.03)
lp_sic_si_3 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb3.nc")
sic_si_3 = get_zeff(lp_sic_si_3, "Si", 0.03, 0.03)
lp_sic_c_2 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb2.nc")
sic_c_2 = get_zeff(lp_sic_c_2, "C", 0.02, 0.02)
lp_sic_si_2 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb2.nc")
sic_si_2 = get_zeff(lp_sic_si_2, "Si", 0.02, 0.02)
lp_sic_c_6 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb6.nc")
sic_c_6 = get_zeff(lp_sic_c_6, "C", 0.06, 0.06)
lp_sic_si_6 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb6.nc")
sic_si_6 = get_zeff(lp_sic_si_6, "Si", 0.06, 0.06)
lp_sic_c_2_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb2_5.nc")
sic_c_2_5 = get_zeff(lp_sic_c_2_5, "C", 0.025, 0.025)
lp_sic_si_2_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb2_5.nc")
sic_si_2_5 = get_zeff(lp_sic_si_2_5, "Si", 0.025, 0.025)
lp_sic_c_3_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb3_5.nc")
sic_c_3_5 = get_zeff(lp_sic_c_3_5, "C", 0.035, 0.035)
lp_sic_si_3_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb3_5.nc")
sic_si_3_5 = get_zeff(lp_sic_si_3_5, "Si", 0.035, 0.035)
lp_sic_c_4_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb4_5.nc")
sic_c_4_5 = get_zeff(lp_sic_c_4_5, "C", 0.045, 0.045)
lp_sic_si_4_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb4_5.nc")
sic_si_4_5 = get_zeff(lp_sic_si_4_5, "Si", 0.045, 0.045)
lp_sic_c_5_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb5_5.nc")
sic_c_5_5 = get_zeff(lp_sic_c_5_5, "C", 0.055, 0.055)
lp_sic_si_5_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb5_5.nc")
sic_si_5_5 = get_zeff(lp_sic_si_5_5, "Si", 0.055, 0.055)

# Runs with a 5 cm gap instead of 10 cm. rlim and rwall are different.
lp_sic_c_4_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb4-5cm-gap.nc")
sic_c_4_5g = get_zeff(lp_sic_c_4_5g, "C", 0.04, 0.04, rlim=0.05, rwall=0.07)
lp_sic_si_4_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb4-5cm-gap.nc")
sic_si_4_5g = get_zeff(lp_sic_si_4_5g, "Si", 0.04, 0.04, rlim=0.05, rwall=0.07)
lp_sic_c_3_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb3-5cm-gap.nc")
sic_c_3_5g = get_zeff(lp_sic_c_3_5g, "C", 0.03, 0.03, rlim=0.05, rwall=0.07)
lp_sic_si_3_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb3-5cm-gap.nc")
sic_si_3_5g = get_zeff(lp_sic_si_4_5g, "Si", 0.03, 0.03, rlim=0.05, rwall=0.07)
lp_sic_c_2_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb2-5cm-gap.nc")
sic_c_2_5g = get_zeff(lp_sic_c_2_5g, "C", 0.02, 0.02, rlim=0.05, rwall=0.07)
lp_sic_si_2_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb2-5cm-gap.nc")
sic_si_2_5g = get_zeff(lp_sic_si_2_5g, "Si", 0.02, 0.02, rlim=0.05, rwall=0.07)
lp_sic_c_5_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb5-5cm-gap.nc")
sic_c_5_5g = get_zeff(lp_sic_c_5_5g, "C", 0.05, 0.05, rlim=0.05, rwall=0.07)
lp_sic_si_5_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb5-5cm-gap.nc")
sic_si_5_5g = get_zeff(lp_sic_si_5_5g, "Si", 0.05, 0.05, rlim=0.05, rwall=0.07)
lp_sic_c_6_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb6-5cm-gap.nc")
sic_c_6_5g = get_zeff(lp_sic_c_6_5g, "C", 0.06, 0.06, rlim=0.05, rwall=0.07)
lp_sic_si_6_5g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb6-5cm-gap.nc")
sic_si_6_5g = get_zeff(lp_sic_si_6_5g, "Si", 0.06, 0.06, rlim=0.05, rwall=0.07)

# Runs with a 2 cm gap instead of 10 cm. rlim and rwall are different.
lp_sic_c_4_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb4-2cm-gap.nc")
sic_c_4_2g = get_zeff(lp_sic_c_4_2g, "C", 0.04, 0.04, rlim=0.02, rwall=0.04)
lp_sic_si_4_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb4-2cm-gap.nc")
sic_si_4_2g = get_zeff(lp_sic_si_4_2g, "Si", 0.04, 0.04, rlim=0.02, rwall=0.04)
lp_sic_c_3_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb3-2cm-gap.nc")
sic_c_3_2g = get_zeff(lp_sic_c_3_2g, "C", 0.03, 0.03, rlim=0.02, rwall=0.04)
lp_sic_si_3_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb3-2cm-gap.nc")
sic_si_3_2g = get_zeff(lp_sic_si_4_2g, "Si", 0.03, 0.03, rlim=0.02, rwall=0.04)
lp_sic_c_2_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb2-2cm-gap.nc")
sic_c_2_2g = get_zeff(lp_sic_c_2_2g, "C", 0.02, 0.02, rlim=0.02, rwall=0.04)
lp_sic_si_2_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb2-2cm-gap.nc")
sic_si_2_2g = get_zeff(lp_sic_si_2_2g, "Si", 0.02, 0.02, rlim=0.02, rwall=0.04)
lp_sic_c_5_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb5-2cm-gap.nc")
sic_c_5_2g = get_zeff(lp_sic_c_5_2g, "C", 0.05, 0.05, rlim=0.02, rwall=0.04)
lp_sic_si_5_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb5-2cm-gap.nc")
sic_si_5_2g = get_zeff(lp_sic_si_5_2g, "Si", 0.05, 0.05, rlim=0.02, rwall=0.04)
lp_sic_c_6_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-c-lamb6-2cm-gap.nc")
sic_c_6_2g = get_zeff(lp_sic_c_6_2g, "C", 0.06, 0.06, rlim=0.02, rwall=0.04)
lp_sic_si_6_2g = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-sic-si-lamb6-2cm-gap.nc")
sic_si_6_2g = get_zeff(lp_sic_si_6_2g, "Si", 0.06, 0.06, rlim=0.02, rwall=0.04)

# The W runs.
# lp_w_6 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb6.nc")
# w_6 = get_zeff(lp_w_6, "W", 0.06, 0.06, sic=False, ftot2_1=8.643e18, ftot2_2=8.643e18, fytot2_1=2.586e17, fytot2_2=2.586e17)

# W Runs.
include_w = False
if include_w:
    lp_w_6 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb6.nc")
    w_6 = get_zeff(lp_w_6, "W", 0.06, 0.06, sic=False)
    lp_w_5_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb5_5.nc")
    w_5_5 = get_zeff(lp_w_5_5, "W", 0.055, 0.055, sic=False)
    lp_w_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb5.nc")
    w_5 = get_zeff(lp_w_5, "W", 0.05, 0.05, sic=False)
    lp_w_4_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb4_5.nc")
    w_4_5 = get_zeff(lp_w_4_5, "W", 0.045, 0.045, sic=False)
    lp_w_4 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb4.nc")
    w_4 = get_zeff(lp_w_4, "W", 0.04, 0.04, sic=False)
    lp_w_3_5 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb3_5.nc")
    w_3_5 = get_zeff(lp_w_3_5, "W", 0.035, 0.035, sic=False)
    lp_w_3 = LimPlots.LimPlots("/Users/zamperini/Documents/d3d_work/lim_runs/1dlim/1dlim-w-lamb3.nc")
    w_3 = get_zeff(lp_w_3, "W", 0.03, 0.03, sic=False)


# Combine the Zeff contributions together for each pair of runs.
r = sic_c_5["rmrs"]
zeff6 = sic_c_6["zeff"] + sic_si_6["zeff"]
zeff5 = sic_c_5["zeff"] + sic_si_5["zeff"]
zeff4 = sic_c_4["zeff"] + sic_si_4["zeff"]
zeff3 = sic_c_3["zeff"] + sic_si_3["zeff"]
zeff2 = sic_c_2["zeff"] + sic_si_2["zeff"]
zeff5_5 = sic_c_5_5["zeff"] + sic_si_5_5["zeff"]
zeff4_5 = sic_c_4_5["zeff"] + sic_si_4_5["zeff"]
zeff3_5 = sic_c_3_5["zeff"] + sic_si_3_5["zeff"]
zeff2_5 = sic_c_2_5["zeff"] + sic_si_2_5["zeff"]

# 5 cm gap cases.
r_5g = sic_c_4_5g["rmrs"]
zeff4_5g = sic_c_4_5g["zeff"] + sic_si_4_5g["zeff"]
zeff3_5g = sic_c_3_5g["zeff"] + sic_si_3_5g["zeff"]
zeff2_5g = sic_c_2_5g["zeff"] + sic_si_2_5g["zeff"]
zeff5_5g = sic_c_5_5g["zeff"] + sic_si_5_5g["zeff"]
zeff6_5g = sic_c_6_5g["zeff"] + sic_si_6_5g["zeff"]

# 2 cm gap cases.
r_2g = sic_c_4_2g["rmrs"]
zeff4_2g = sic_c_4_2g["zeff"] + sic_si_4_2g["zeff"]
zeff3_2g = sic_c_3_2g["zeff"] + sic_si_3_2g["zeff"]
zeff2_2g = sic_c_2_2g["zeff"] + sic_si_2_2g["zeff"]
zeff5_2g = sic_c_5_2g["zeff"] + sic_si_5_2g["zeff"]
zeff6_2g = sic_c_6_2g["zeff"] + sic_si_6_2g["zeff"]

# Graphite comparisons.
zeff6_gph = sic_c_6["zeff_gph"]
zeff5_gph = sic_c_5["zeff_gph"]
zeff4_gph = sic_c_4["zeff_gph"]
zeff3_gph = sic_c_3["zeff_gph"]
zeff2_gph = sic_c_2["zeff_gph"]
zeff5_5_gph = sic_c_5_5["zeff_gph"]
zeff4_5_gph = sic_c_4_5["zeff_gph"]
zeff3_5_gph = sic_c_3_5["zeff_gph"]
zeff2_5_gph = sic_c_2_5["zeff_gph"]

# 5 cm gap cases.
zeff4_5g_gph = sic_c_4_5g["zeff_gph"]
zeff3_5g_gph = sic_c_3_5g["zeff_gph"]
zeff2_5g_gph = sic_c_2_5g["zeff_gph"]
zeff5_5g_gph = sic_c_5_5g["zeff_gph"]
zeff6_5g_gph = sic_c_6_5g["zeff_gph"]

# 2 cm gap cases.
zeff4_2g_gph = sic_c_4_2g["zeff_gph"]
zeff3_2g_gph = sic_c_3_2g["zeff_gph"]
zeff2_2g_gph = sic_c_2_2g["zeff_gph"]
zeff5_2g_gph = sic_c_5_2g["zeff_gph"]
zeff6_2g_gph = sic_c_6_2g["zeff_gph"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

ax1.plot(r * 100, zeff6, color="C1")
ax1.plot(r * 100, zeff5, color="C2")
ax1.plot(r * 100, zeff4, color="C3")
ax1.plot(r * 100, zeff3, color="C4")
ax1.plot(r * 100, zeff2, color="C5")
if include_w:
    ax1.plot(r * 100, w_6["zeff"], color="C1", linestyle="--")
    ax1.plot(r * 100, w_3["zeff"], color="C4", linestyle="--")
ax1.set_xlabel(r"$\mathdefault{R-R_{sep}}$ (cm)")
ax1.set_ylabel(r"$\mathdefault{Z_{eff}}$ Contribution", fontsize=14)

# Values at the separatrix.
ax2.scatter(6, zeff6[-1], color="tab:red")
ax2.scatter(5, zeff5[-1], color="tab:red")
ax2.scatter(4, zeff4[-1], color="tab:red")
ax2.scatter(3, zeff3[-1], color="tab:red")
ax2.scatter(2, zeff2[-1], color="tab:red")
ax2.scatter(5.5, zeff5_5[-1], color="tab:red")
ax2.scatter(4.5, zeff4_5[-1], color="tab:red")
ax2.scatter(3.5, zeff3_5[-1], color="tab:red")
ax2.scatter(2.5, zeff2_5[-1], color="tab:red")
ax2.scatter(6, zeff6_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(5, zeff5_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(4, zeff4_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(3, zeff3_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(2, zeff2_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(5.5, zeff5_5_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(4.5, zeff4_5_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(3.5, zeff3_5_gph[-1], edgecolors="tab:red", color="w")
ax2.scatter(2.5, zeff2_5_gph[-1], edgecolors="tab:red", color="w")
if include_w:
    ax2.scatter(6, w_6["zeff"][-1], color="tab:red", marker="^", edgecolors="k")
    ax2.scatter(5, w_5["zeff"][-1], color="tab:red", marker="^", edgecolors="k")
    ax2.scatter(4, w_4["zeff"][-1], color="tab:red", marker="^", edgecolors="k")
    ax2.scatter(3, w_3["zeff"][-1], color="tab:red", marker="^", edgecolors="k")
ax2.set_xlabel(r"$\lambda$ = $\lambda_n$ = $\lambda_{T}$ (cm)", fontsize=14)
ax2.set_ylabel(r"$\mathdefault{Z_{eff}}$ Contribution @ Separatrix", fontsize=14)

fig.tight_layout()
fig.show()

# Easier to just do a new plot with lines of Zeff at the separatrix of each gap case.
lw = 3
fig, ax1 = plt.subplots(figsize=(5, 4))

# 2 cm gaps.
lambdas = [2.0, 3.0, 4.0, 5.0, 6.0]
all_sic_zeffs = [zeff2_2g[-1], zeff3_2g[-1], zeff4_2g[-1], zeff5_2g[-1], zeff6_2g[-1]]
all_gph_zeffs = [zeff2_2g_gph[-1], zeff3_2g_gph[-1], zeff4_2g_gph[-1], zeff5_2g_gph[-1], zeff6_2g_gph[-1]]
ax1.plot(lambdas, all_sic_zeffs, color="tab:green", label="2 cm", lw=lw)
ax1.plot(lambdas, all_gph_zeffs, color="tab:green", linestyle="--", lw=lw)

# 5 cm gaps.
lambdas = [2.0, 3.0, 4.0, 5.0, 6.0]
all_sic_zeffs = [zeff2_5g[-1], zeff3_5g[-1], zeff4_5g[-1], zeff5_5g[-1], zeff6_5g[-1]]
all_gph_zeffs = [zeff2_5g_gph[-1], zeff3_5g_gph[-1], zeff4_5g_gph[-1], zeff5_5g_gph[-1], zeff6_5g_gph[-1]]
ax1.plot(lambdas, all_sic_zeffs, color="tab:red", label="5 cm", lw=lw)
ax1.plot(lambdas, all_gph_zeffs, color="tab:red", linestyle="--", lw=lw)

# 10 cm gaps.
lambdas = np.arange(2.0, 6.5, 0.5)
all_sic_zeffs = [zeff2[-1], zeff2_5[-1], zeff3[-1], zeff3_5[-1], zeff4[-1], zeff4_5[-1], zeff5[-1], zeff5_5[-1],
                 zeff6[-1]]
all_gph_zeffs = [zeff2_gph[-1], zeff2_5_gph[-1], zeff3_gph[-1], zeff3_5_gph[-1], zeff4_gph[-1], zeff4_5_gph[-1],
                 zeff5_gph[-1], zeff5_5_gph[-1], zeff6_gph[-1]]
ax1.plot(lambdas, all_sic_zeffs, color="tab:purple", label="10 cm", lw=lw)
ax1.plot(lambdas, all_gph_zeffs, color="tab:purple", linestyle="--", lw=lw)

ax1.set_xlabel(r"$\lambda$ = $\lambda_n$ = $\lambda_{T}$ (cm)", fontsize=16)
ax1.set_ylabel(r"$\mathdefault{Z_{eff}}$ Contribution at Sep.", fontsize=16)
ax1.legend(fontsize=14)
ax1.set_yticks([0, 0.5, 1, 1.5])
ax1.set_ylim([0, None])
ax1.tick_params(axis="both", labelsize=14)

fig.tight_layout()
fig.show()
