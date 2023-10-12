# This script is based off divimp_zeff_sic_range.py. Except now we have multiple cases that used varying backgrounds,
# and we want to compare how the relative improvement in Zeff is when the target Te is increased to more intense
# # conditions.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

def calc_zeffs(grap_path, si_path, c_path):

    # 009: Normal case
    # 012: Background target multipliers set to 2.0 (mult6).
    grap = oedge_plots.OedgePlots(grap_path)

    # 009: fC = 2%, fSi = 0.2%
    # 012: Same but with target multipliers for background set to 2.0 (mult6).
    si_max = oedge_plots.OedgePlots(si_path)
    c_max = oedge_plots.OedgePlots(c_path)

    # 011: fC = 0.2%, fSi = 0.02%
    # 013: Same but with target multipliers for background set to 2.0 (mult6).
    #si_min = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-013.nc")
    #c_min = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-013.nc")

    dsum1 = 0.0
    dsum2 = 0.0
    dsum1_c_max = 0.0
    dsum2_c_max = 0.0
    dsum1_si_max = 0.0
    dsum2_si_max = 0.0
    dsum1_c_min = 0.0
    dsum2_c_min = 0.0
    dsum1_si_min = 0.0
    dsum2_si_min = 0.0
    dsum1_c6 = 0.0
    dsum2_c6 = 0.0

    # Calculate Zeff using the way DIVIMP does it. I admit I don't fully understand it, but obvi DIVIMP is trustworthy.
    for charge in range(1, 30):
        print("Charge: {}".format(charge))
        if charge < 7:
            #nz = np.array(c_min.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
            #dsum1_c_min += nz * charge
            #dsum2_c_min += nz * charge ** 2

            nz = np.array(c_max.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
            dsum1_c_max += nz * charge
            dsum2_c_max += nz * charge ** 2

            nz = np.array(grap.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
            dsum1 += nz * charge
            dsum2 += nz * charge ** 2

            # The CER Zeff calculation is assuming C6+ is the only impurity charge
            # state that exists, so that's what we want to compare to.
            if charge == 6:
                dsum1_c6 += nz * charge
                dsum2_c6 += nz * charge ** 2

        # Si only goes up to 14.
        if charge < 15:
            #nz = np.array(si_min.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
            #dsum1_si_min += nz * charge
            #dsum2_si_min += nz * charge ** 2

            nz = np.array(si_max.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
            dsum1_si_max += nz * charge
            dsum2_si_max += nz * charge ** 2

    # Needed for below.
    ne_probe = grap.fake_probe(1.9, 2.5, 0.0, 0.0, "ne")
    divimp_r = np.array(ne_probe["r"])
    ne = np.array(ne_probe["ne"])
    zeff_sic = np.zeros(ne.shape)
    zeff_gph = np.zeros(ne.shape)

    # Mimicking calculation from div.f.
    zeffs1 = dsum1
    zeffs2 = 1.0 * ne - zeffs1
    zeff_gph = np.zeros(zeffs1.shape)
    for i in range(0, len(zeff_gph)):
        if zeffs2[i] > 0:
            zeff_gph[i] = (1.0 * zeffs2[i] + dsum2[i]) / (1.0 * ne[i])
        else:
            zeff_gph[i] = dsum2[i] / dsum1[i]

    # Again for SiC, min and max.
    # dsum1_sic = dsum1_si_min + dsum1_c_min
    # dsum2_sic = dsum2_si_min + dsum2_c_min
    # zeffs1 = dsum1_sic
    # zeffs2 = 1.0 * ne - zeffs1
    # zeff_sic_min = np.zeros(len(zeffs1))
    # for i in range(0, len(zeff_sic)):
    #     if zeffs2[i] > 0:
    #         zeff_sic_min[i] = (1.0 * zeffs2[i] + dsum2_sic[i]) / (1.0 * ne[i])
    #     else:
    #         zeff_sic_min[i] = dsum2_sic[i] / dsum1_sic[i]

    dsum1_sic = dsum1_si_max + dsum1_c_max
    dsum2_sic = dsum2_si_max + dsum2_c_max
    zeffs1 = dsum1_sic
    zeffs2 = 1.0 * ne - zeffs1
    zeff_sic_max = np.zeros(len(zeffs1))
    for i in range(0, len(zeff_sic)):
        if zeffs2[i] > 0:
            zeff_sic_max[i] = (1.0 * zeffs2[i] + dsum2_sic[i]) / (1.0 * ne[i])
        else:
            zeff_sic_max[i] = dsum2_sic[i] / dsum1_sic[i]

    # One more time for the charge state being viewed.
    zeffs1 = dsum1_c6
    zeffs2 = 1.0 * ne - zeffs1
    zeff_gph_c6 = np.zeros(len(zeffs1))
    for i in range(0, len(zeff_gph_c6)):
        if zeffs2[i] > 0:
            zeff_gph_c6[i] = (1.0 * zeffs2[i] + dsum2_c6[i]) / (1.0 * ne[i])
        else:
            zeff_gph_c6[i] = dsum2_c6[i] / dsum1_c6[i]

    return {"divimp_r":divimp_r, "zeff_sic_max":zeff_sic_max, "zeff_gph":zeff_gph}

#Load the Zeffs for each scenario.
# The OG, 009 case with the constrained background. Targ Te mults = 0.50.
zeff_og = calc_zeffs("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-009.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-009.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-009.nc")
# The case with mult1, Targ Te mults = 0.75.
zeff_mult1 = calc_zeffs("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-014.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-014.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-014.nc")
# The case with mult1, Targ Te mults = 1.00.
zeff_mult2 = calc_zeffs("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-015.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-015.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-015.nc")
# The case with mult1, Targ Te mults = 1.25.
zeff_mult3 = calc_zeffs("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-016.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-016.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-016.nc")
# The case with mult1, Targ Te mults = 1.50.
zeff_mult4 = calc_zeffs("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-017.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-017.nc",
                     "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-017.nc")


cmap = plt.get_cmap('inferno')
colors = cmap(np.linspace(0, 0.9, 5))
r = zeff_og["divimp_r"]
divmask = r > 2.21
y0 = savgol_filter(zeff_og["zeff_sic_max"][divmask], 5, 2) / savgol_filter(zeff_og["zeff_gph"][divmask], 5, 2)
y1 = savgol_filter(zeff_mult1["zeff_sic_max"][divmask], 5, 2) / savgol_filter(zeff_mult1["zeff_gph"][divmask], 5, 2)
y2 = savgol_filter(zeff_mult2["zeff_sic_max"][divmask], 5, 2) / savgol_filter(zeff_mult2["zeff_gph"][divmask], 5, 2)
y3 = savgol_filter(zeff_mult3["zeff_sic_max"][divmask], 5, 2) / savgol_filter(zeff_mult3["zeff_gph"][divmask], 5, 2)
y4 = savgol_filter(zeff_mult4["zeff_sic_max"][divmask], 5, 2) / savgol_filter(zeff_mult4["zeff_gph"][divmask], 5, 2)
mask = r > 2.26
lw = 3
fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.axvline(2.26, color="k", linestyle="--")
ax1.axhline(1.00, color="k")
ax1.plot(r[divmask], y0, lw=lw, label="Base", color=colors[0])
ax1.plot(r[divmask], y1, lw=lw, label="1.5x", color=colors[1])
ax1.plot(r[divmask], y2, lw=lw, label="2.0x", color=colors[2])
ax1.plot(r[divmask], y3, lw=lw, label="2.25x", color=colors[3])
ax1.plot(r[divmask], y4, lw=lw, label="2.5x", color=colors[4])
ax1.set_xlabel("R (m)", fontsize=12)
ax1.set_ylabel(r"$\mathdefault{Z_{eff}^{SiC}\ /\ Z_{eff}^{Graphite}}$", fontsize=12)
ax1.legend()
ax1.tick_params(axis="both", labelsize=10)
fig.tight_layout()
fig.show()