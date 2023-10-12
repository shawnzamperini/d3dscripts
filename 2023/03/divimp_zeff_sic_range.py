# Script to plot the Zeff comparison from SiC as compared ot graphite for a range of the fractions for the Si and C
# # sputtering values in the SiC mixed material model.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from gadata import gadata
import MDSplus

# 009: Normal case
# 012: Background target multipliers set to 2.0 (mult6).
grap = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-012.nc")

# 009: fC = 2%, fSi = 0.2%
# 012: Same but with target multipliers for background set to 2.0 (mult6).
si_max = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-012.nc")
c_max = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-012.nc")

# 011: fC = 0.2%, fSi = 0.02%
# 013: Same but with target multipliers for background set to 2.0 (mult6).
si_min = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-013.nc")
c_min = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-013.nc")

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
        nz = np.array(c_min.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        dsum1_c_min += nz * charge
        dsum2_c_min += nz * charge ** 2

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
        nz = np.array(si_min.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        dsum1_si_min += nz * charge
        dsum2_si_min += nz * charge ** 2

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
dsum1_sic = dsum1_si_min + dsum1_c_min
dsum2_sic = dsum2_si_min + dsum2_c_min
zeffs1 = dsum1_sic
zeffs2 = 1.0 * ne - zeffs1
zeff_sic_min = np.zeros(len(zeffs1))
for i in range(0, len(zeff_sic)):
    if zeffs2[i] > 0:
        zeff_sic_min[i] = (1.0 * zeffs2[i] + dsum2_sic[i]) / (1.0 * ne[i])
    else:
        zeff_sic_min[i] = dsum2_sic[i] / dsum1_sic[i]

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

# Need the densities for C6+ and all carbon to compare.
divimp_nz_all = np.array(grap.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"])
divimp_nz_c6 = np.array(grap.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=6)["nz"])

tchords = ["T{}".format(i) for i in range(1, 49)]
vchords = ["V{}".format(i) for i in range(1, 33)]
chords = np.append(tchords, vchords)

conn = MDSplus.Connection("atlas.gat.com")
shot = 167196
tmin = 2500
tmax = 5000
cer = {}
for chord in chords:
    print("Chord: {}".format(chord))
    signal_r = "CERAR{}".format(chord)
    signal_z = "CERAZ{}".format(chord)
    signal_zeff = "CERAZEFF{}".format(chord)
    signal_nz = "CERANZ{}".format(chord)
    gaobj_r = gadata(signal_r, shot, connection=conn)
    gaobj_z = gadata(signal_z, shot, connection=conn)
    gaobj_zeff = gadata(signal_zeff, shot, connection=conn)
    gaobj_nz = gadata(signal_nz, shot, connection=conn)
    chord_r = float(gaobj_r.zdata[0])
    chord_z = float(gaobj_z.zdata[0])

    # Return average Zeff value.
    mask = np.logical_and(gaobj_zeff.xdata > tmin, gaobj_zeff.xdata < tmax)
    chord_zeff = float(gaobj_zeff.zdata[mask].mean())
    chord_nz = float(gaobj_nz.zdata[mask].mean())
    chord_zeff_std = float(gaobj_zeff.zdata[mask].std())
    chord_nz_std = float(gaobj_nz.zdata[mask].std())

    # If zero, then no data for this chord.
    if chord_zeff != 0:
        cer[chord] = {"r": chord_r, "z": chord_z, "zeff": chord_zeff,
                      "nz": chord_nz, "zeff_err": chord_zeff_std, "nz_err": chord_nz_std}

# Distill into simple arrays.
tchords_r = [cer[c]["r"] for c in cer.keys() if c[0] == "T"]
tchords_zeff = [cer[c]["zeff"] for c in cer.keys() if c[0] == "T"]
tchords_nz = [cer[c]["nz"] for c in cer.keys() if c[0] == "T"]
tchords_zeff_err = [cer[c]["zeff_err"] for c in cer.keys() if c[0] == "T"]
tchords_nz_err = [cer[c]["nz_err"] for c in cer.keys() if c[0] == "T"]
vchords_r = [cer[c]["r"] for c in cer.keys() if c[0] == "V"]
vchords_zeff = [cer[c]["zeff"] for c in cer.keys() if c[0] == "V"]
vchords_nz = [cer[c]["nz"] for c in cer.keys() if c[0] == "V"]
vchords_zeff_err = [cer[c]["zeff_err"] for c in cer.keys() if c[0] == "V"]
vchords_nz_err = [cer[c]["nz_err"] for c in cer.keys() if c[0] == "V"]

divmask = divimp_r > 2.21
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.axvline(2.26, color="k", linestyle="--")
ax1.plot(divimp_r[divmask], zeff_gph[divmask], label="Graphite", color="tab:red", lw=3)
ax1.fill_between(divimp_r[divmask], zeff_sic_min[divmask], zeff_sic_max[divmask], label="SiC", alpha=0.4,
                 color="tab:green")

ax1.legend(fontsize=14)
ax1.set_xlabel("R (m)", fontsize=16)
ax1.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=16)
ax1.tick_params(axis="both", labelsize=14)
ax1.set_ylim([1.0, 2.0])
fig.tight_layout()
fig.show()

# Again but with all three plots together for 4 pager SiC paper.
fontsize = 12
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 3), sharex=True)

ax1.axvline(2.26, color="k", linestyle="--")
ax1.scatter(tchords_r, tchords_nz, c="tab:red", marker="o", label="CER", edgecolors="k")
# ax1.scatter(vchords_r, vchords_nz, c="tab:red", marker="*", label="CER (V)", edgecolors="k")
# ax1.plot(divimp_r, divimp_nz_all, color="tab:red", label="DIVIMP-Graphite")
ax1.plot(divimp_r[divmask], divimp_nz_c6[divmask], color="k", linestyle="-", lw=3)
ax1.plot(divimp_r[divmask], divimp_nz_c6[divmask], color="tab:red", linestyle="-", label="DIVIMP", lw=2)
ax1.legend(fontsize=fontsize-2, framealpha=1.0)
ax1.set_xlabel("R (m)", fontsize=fontsize-2)
ax1.set_ylabel(r"C6+ Density ($\mathdefault{m^{-3}}$)", fontsize=fontsize)
ax1.set_ylim(0, 0.8e18)
ax1.set_title("C6+ Only")
ax1.text(0.02, 0.92, "a)", transform=ax1.transAxes)

ax2.axvline(2.26, color="k", linestyle="--")
ax2.scatter(tchords_r, tchords_zeff, c="tab:red", marker="o", label="CER", edgecolors="k")
# ax2.scatter(vchords_r, vchords_zeff, c="tab:red", marker="*", label="CER (V)", edgecolors="k")
# ax2.plot(divimp_r, zeff_gph, color="tab:red", label="Graphite")
ax2.plot(divimp_r[divmask], zeff_gph_c6[divmask], color="k", linestyle="-", lw=3)
ax2.plot(divimp_r[divmask], zeff_gph_c6[divmask], color="tab:red", label="DIVIMP", linestyle="-", lw=2)
# ax2.plot(divimp_r, zeff_gph_code, color="tab:red", linestyle="--", label="Code-calculated")
# ax2.plot(divimp_r, zeff_sic, color="tab:purple", label="SiC")
# ax2.plot(divimp_r, zeff_w, color="tab:green", label="W")
ax2.legend(fontsize=fontsize-2, framealpha=1.0)
ax2.set_ylim([1, 2.5])
ax2.set_xlim([2.2, 2.4])
ax2.set_xlabel("R (m)", fontsize=fontsize)
ax2.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=fontsize)
ax2.set_title("C6+ Only", fontsize=fontsize)
ax2.set_ylim([1, 2.25])
ax2.set_yticks([1.0, 1.5, 2.0])
ax2.text(0.02, 0.92, "b)", transform=ax2.transAxes)

ax3.axvline(2.26, color="k", linestyle="--")
ax3.plot(divimp_r[divmask], zeff_gph[divmask], label="Graphite", color="tab:red", lw=3)
ax3.fill_between(divimp_r[divmask], zeff_sic_min[divmask], zeff_sic_max[divmask], label="SiC", alpha=0.4,
                 color="tab:green")

ax3.legend(fontsize=fontsize-2, framealpha=1.0)
ax3.set_xlabel("R (m)", fontsize=fontsize)
ax3.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=fontsize)
# ax3.tick_params(axis="both", labelsize=14)
ax3.set_ylim([1.0, 2.25])
ax3.set_title("All Charge States")
ax3.set_yticks([1.0, 1.5, 2.0])
ax3.text(0.02, 0.92, "c)", transform=ax3.transAxes)

fig.tight_layout()
fig.show()

# Rough estimate of the fuel dilution.
grap_dil = (ne_probe["ne"] - dsum1) / ne_probe["ne"]
sic_dil_min = (ne_probe["ne"] - dsum1_c_min - dsum1_si_min) / ne_probe["ne"]
sic_dil_max = (ne_probe["ne"] - dsum1_c_max - dsum1_si_max) / ne_probe["ne"]
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.axvline(2.26, color="k", linestyle="--")
ax1.plot(divimp_r[divmask], grap_dil[divmask], label="Graphite", color="tab:red", lw=3)
ax1.fill_between(divimp_r[divmask], sic_dil_min[divmask], sic_dil_max[divmask], label="SiC", alpha=0.4,
                 color="tab:green")

ax1.legend(fontsize=14)
ax1.set_xlabel("R (m)", fontsize=16)
ax1.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=16)
ax1.tick_params(axis="both", labelsize=14)
#ax1.set_ylim([1.0, 2.0])
fig.tight_layout()
fig.show()