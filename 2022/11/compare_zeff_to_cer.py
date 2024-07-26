# Script to load CER data and compare the DIVIMP simulation for 167196 to it.
import matplotlib.pyplot as plt
import numpy as np
from gadata import gadata
import oedge_plots
import MDSplus

# 005: The original baseline scenario with diffusion (Dperp = 1.0 m2/s).
# 006: Used a blobby model with 167195 blob data.
# 008: Dperp = 0.3 m2/s
# 009: Core Dperp = 0.3, SOL Dperp = 1.0.
cpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-c-009.nc"
spath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-allsic-wall-si-009.nc"
gpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allgrap-009.nc"
# wpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allw-2pc-006.nc"
wpath = "/Users/zamperini/Documents/d3d_work/divimp_files/sic_wall/d3d-sic-wall-allw-009.nc"
cop = oedge_plots.OedgePlots(cpath)
sop = oedge_plots.OedgePlots(spath)
gop = oedge_plots.OedgePlots(gpath)
wop = oedge_plots.OedgePlots(wpath)

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

# Calculate Zeff from DIVIMP at the CER location.
ne_probe = cop.fake_probe(1.9, 2.5, 0.0, 0.0, "ne")
divimp_r = np.array(ne_probe["r"])
ne = np.array(ne_probe["ne"])
zeff_sic = np.zeros(ne.shape)
zeff_gph = np.zeros(ne.shape)
cop_absfac = float(cop.absfac)
gop_absfac = float(gop.absfac)
sop_absfac = float(sop.absfac)
wop_absfac = float(wop.absfac)
dsum1 = np.zeros(ne.shape)
dsum2 = np.zeros(ne.shape)
dsum1_c6 = np.zeros(ne.shape)
dsum2_c6 = np.zeros(ne.shape)
dsum1_si = np.zeros(ne.shape)
dsum2_si = np.zeros(ne.shape)
dsum1_c = np.zeros(ne.shape)
dsum2_c = np.zeros(ne.shape)
dsum1_w = np.zeros(ne.shape)
dsum2_w = np.zeros(ne.shape)
for charge in range(1, 30):
    print("Charge: {}".format(charge))
    if charge < 7:
        nz = np.array(cop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        # zeff_sic += charge**2 * nz / ne
        dsum1_c += nz * charge
        dsum2_c += nz * charge ** 2

        nz = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        # zeff_gph += charge**2 * nz / ne

        # Let's just mimic the calculation in div.f.
        dsum1 += nz * charge
        dsum2 += nz * charge ** 2

        # The CER Zeff calculation is assuming C6+ is the only impurity charge
        # state that exists, so that's what we want to compare to.
        if charge == 6:
            dsum1_c6 += nz * charge
            dsum2_c6 += nz * charge ** 2

    # Si only goes up to 14...
    if charge < 15:
        nz = np.array(sop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
        print()
        # zeff_sic += charge**2 * nz / ne
        dsum1_si += nz * charge
        dsum2_si += nz * charge ** 2

    # And then W.
    nz = np.array(wop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"])
    dsum1_w += nz * charge
    dsum2_w += nz * charge ** 2

# Mimicking calculation from div.f.
prompt_redep_frac = 0.0
zeffs1 = dsum1 * (1 - prompt_redep_frac)
zeffs2 = 1.0 * ne - zeffs1
zeff_gph = np.zeros(zeffs1.shape)
for i in range(0, len(zeff_gph)):
    if zeffs2[i] > 0:
        zeff_gph[i] = (1.0 * zeffs2[i] + dsum2[i] * (1 - prompt_redep_frac)) / (1.0 * ne[i])
    else:
        zeff_gph[i] = dsum2[i] * (1 - prompt_redep_frac) / (dsum1[i] * (1 - prompt_redep_frac))

# Again for SiC.
dsum1_sic = dsum1_si + dsum1_c
dsum2_sic = dsum2_si + dsum2_c
zeffs1 = dsum1_sic * (1 - prompt_redep_frac)
zeffs2 = 1.0 * ne - zeffs1
zeff_sic = np.zeros(len(zeffs1))
for i in range(0, len(zeff_sic)):
    if zeffs2[i] > 0:
        zeff_sic[i] = (1.0 * zeffs2[i] + dsum2_sic[i] * (1 - prompt_redep_frac)) / (1.0 * ne[i])
    else:
        zeff_sic[i] = dsum2_sic[i] * (1 - prompt_redep_frac) / (dsum1_sic[i] * (1 - prompt_redep_frac))

# One more time for the charge state being viewed.
zeffs1 = dsum1_c6 * (1 - prompt_redep_frac)
zeffs2 = 1.0 * ne - zeffs1
zeff_gph_c6 = np.zeros(len(zeffs1))
for i in range(0, len(zeff_gph_c6)):
    if zeffs2[i] > 0:
        zeff_gph_c6[i] = (1.0 * zeffs2[i] + dsum2_c6[i] * (1 - prompt_redep_frac)) / (1.0 * ne[i])
    else:
        zeff_gph_c6[i] = dsum2_c6[i] * (1 - prompt_redep_frac) / (dsum1_c6[i] * (1 - prompt_redep_frac))

zeffs1 = dsum1_w * (1 - prompt_redep_frac)
zeffs2 = 1.0 * ne - zeffs1
zeff_w = np.zeros(zeffs1.shape)
for i in range(0, len(zeff_w)):
    if zeffs2[i] > 0:
        zeff_w[i] = (1.0 * zeffs2[i] + dsum2_w[i] * (1 - prompt_redep_frac)) / (1.0 * ne[i])
    else:
        zeff_w[i] = dsum2_w[i] * (1 - prompt_redep_frac) / (dsum1_w[i] * (1 - prompt_redep_frac))

# Need the densities for C6+ and all carbon to compare.
divimp_nz_all = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"])
divimp_nz_c6 = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=6)["nz"])

# Apply any prompt redeposition fudge factor, and then plus 1 for deuterium.
# zeff_sic_only = zeff_sic * (1-prompt_redep_frac)
# zeff_gph_only = zeff_gph * (1-prompt_redep_frac)
# zeff_sic_final = zeff_sic_only  # Took out the +1, don't think it's right.
# zeff_gph_final = zeff_gph_only

# Comparison to the code calculated value.
# zeff_gph_code = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "zeff")["zeff"])

# The impurity radiated power. Maybe not super useful since W has the additional thing that it radiates Bremstrahlung
# as well, which I don't think is included here.
powl_s = np.array(sop.fake_probe(1.9, 2.5, 0.0, 0.0, "prad", charge="all")["prad"])
powl_c = np.array(cop.fake_probe(1.9, 2.5, 0.0, 0.0, "prad", charge="all")["prad"])
powl_g = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "prad", charge="all")["prad"])
powl_w = np.array(wop.fake_probe(1.9, 2.5, 0.0, 0.0, "prad", charge="all")["prad"])
ns = np.array(sop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"])
nc = np.array(cop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"])
ng = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"])
nw = np.array(wop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"])
vol = np.array(sop.fake_probe(1.9, 2.5, 0.0, 0.0, "vol")["vol"])
prad_s = powl_s * ns * vol
prad_c = powl_c * nc * vol
prad_g = powl_g * ng * vol
prad_w = powl_w * nw * vol

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

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 4), sharex=True)

ax2.axvline(2.26, color="k", linestyle="--")
ax2.scatter(tchords_r, tchords_zeff, c="tab:red", marker="o", label="T")
ax2.scatter(vchords_r, vchords_zeff, c="tab:red", marker="*", label="V")
ax2.plot(divimp_r, zeff_gph, color="tab:red", label="Graphite")
ax2.plot(divimp_r, zeff_gph_c6, color="tab:red", label="Graphite (C6+)", linestyle="--")
# ax2.plot(divimp_r, zeff_gph_code, color="tab:red", linestyle="--", label="Code-calculated")
ax2.plot(divimp_r, zeff_sic, color="tab:purple", label="SiC")
ax2.plot(divimp_r, zeff_w, color="tab:green", label="W")
ax2.legend()
ax2.set_ylim([1, 3])
ax2.set_xlim([2.1, 2.4])
ax2.set_xlabel("R (m)")
ax2.set_ylabel("Zeff")

ax1.axvline(2.26, color="k", linestyle="--")
ax1.scatter(tchords_r, tchords_nz, c="tab:red", marker="o", label="T")
ax1.scatter(vchords_r, vchords_nz, c="tab:red", marker="*", label="V")
ax1.plot(divimp_r, divimp_nz_all, color="tab:red", label="DIVIMP-Graphite")
ax1.plot(divimp_r, divimp_nz_c6, color="tab:red", linestyle="--", label="DIVIMP-Graphite (C6+)")
ax1.legend()
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Carbon Density (m-3)")
ax1.set_ylim(0, 1.5e18)

ax3.axvline(2.26, color="k", linestyle="--")
ax3.plot(divimp_r, prad_g, color="tab:red")
ax3.plot(divimp_r, prad_s + prad_c, color="tab:purple")
ax3.plot(divimp_r, prad_w, color="tab:green")
ax3.set_xlabel("R (m)")
ax3.set_ylabel("Imp. Prad (W/m3)")
ax3.set_yscale("log")
ax3.set_ylim([1e-5, 10])

fig.tight_layout()
fig.show()

# Plot of just the graphite results to show the simulation results should be trustworthy.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)

ax1.axvline(2.26, color="k", linestyle="--")
ax1.scatter(tchords_r, tchords_nz, c="tab:red", marker="o", label="CER", edgecolors="k")
# ax1.scatter(vchords_r, vchords_nz, c="tab:red", marker="*", label="CER (V)", edgecolors="k")
ax1.plot(divimp_r, divimp_nz_all, color="tab:red", label="DIVIMP-Graphite")
ax1.plot(divimp_r[divmask], divimp_nz_c6[divmask], color="k", linestyle="-", lw=3)
ax1.plot(divimp_r[divmask], divimp_nz_c6[divmask], color="tab:red", linestyle="--", label="DIVIMP (C6+)", lw=2)
ax1.legend(fontsize=12)
ax1.set_xlabel("R (m)", fontsize=14)
ax1.set_ylabel(r"Carbon Density ($\mathdefault{m^{-3}}$)", fontsize=14)
ax1.set_ylim(0, 1.0e18)

ax2.axvline(2.26, color="k", linestyle="--")
ax2.scatter(tchords_r, tchords_zeff, c="tab:red", marker="o", label="CER", edgecolors="k")
# ax2.scatter(vchords_r, vchords_zeff, c="tab:red", marker="*", label="CER (V)", edgecolors="k")
ax2.plot(divimp_r, zeff_gph, color="tab:red", label="Graphite")
ax2.plot(divimp_r[divmask], zeff_gph_c6[divmask], color="k", linestyle="-", lw=3)
ax2.plot(divimp_r[divmask], zeff_gph_c6[divmask], color="tab:red", label="DIVIMP (C6+)", linestyle="--", lw=2)
# ax2.plot(divimp_r, zeff_gph_code, color="tab:red", linestyle="--", label="Code-calculated")
# ax2.plot(divimp_r, zeff_sic, color="tab:purple", label="SiC")
# ax2.plot(divimp_r, zeff_w, color="tab:green", label="W")
ax2.legend(fontsize=13)
ax2.set_ylim([1, 2.5])
ax2.set_xlim([2.1, 2.4])
ax2.set_xlabel("R (m)", fontsize=14)
ax2.set_ylabel(r"$\mathdefault{Z_{eff}}$", fontsize=14)

fig.tight_layout()
fig.show()