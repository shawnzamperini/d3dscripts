# Script to compare Zeff and nz for C6+ as measured by CER for a 190423 DIVIMP
# carbon run.
import matplotlib.pyplot as plt
import numpy as np
from gadata import gadata
import oedge_plots
import MDSplus


gpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-carbon-002.nc"
gop = oedge_plots.OedgePlots(gpath)

tchords = ["T{}".format(i) for i in range(1, 49)]
vchords = ["V{}".format(i) for i in range(1, 33)]
chords = np.append(tchords, vchords)

conn = MDSplus.Connection("atlas.gat.com")
shot = 190423
tmin = 2800
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
        cer[chord] = {"r":chord_r, "z": chord_z, "zeff":chord_zeff,
            "nz":chord_nz, "zeff_err":chord_zeff_std, "nz_err":chord_nz_std}

# Calculate Zeff from DIVIMP at the CER location.
ne_probe = gop.fake_probe(1.9, 2.5, 0.0, 0.0, "ne")
divimp_r = np.array(ne_probe["r"])
ne = np.array(ne_probe["ne"])
gop_absfac = float(gop.absfac)
dsum1 = 0.0
dsum2 = 0.0
dsum1_c6 = 0.0
dsum2_c6 = 0.0
for charge in range(1, 7):
    print("Charge: {}".format(charge))

    nz = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=charge)["nz"]) * gop_absfac
    dsum1 += nz * charge
    dsum2 += nz * charge**2

    # The CER Zeff calculation is assuming C6+ is the only impurity charge
    # state that exists, so that's what we want to compare to.
    if charge == 6:
        dsum1_c6 += nz * charge
        dsum2_c6 += nz * charge**2

# This modifier is a fudge factor to account for all the carbon DIVIMP is unable
# to track, e.g., from the inner wall outside of the grid.
nz_mod = 750

# Mimicing calculation from div.f.
zeffs1 = dsum1 * nz_mod
zeffs2 = 1.0 * ne - zeffs1
zeff_gph = np.zeros(len(zeffs1))
for i in range(0, len(zeff_gph)):
    if zeffs2[i] > 0:
        zeff_gph[i] = (1.0 * zeffs2[i] + dsum2[i] * nz_mod) / (1.0 * ne[i])
    else:
        zeff_gph[i] = dsum2[i] / dsum1[i]

# One more time for the charge state being viewed.
zeffs1 = dsum1_c6 * nz_mod
zeffs2 = 1.0 * ne - zeffs1
zeff_gph_c6 = np.zeros(len(zeffs1))
for i in range(0, len(zeff_gph_c6)):
    if zeffs2[i] > 0:
        zeff_gph_c6[i] = (1.0 * zeffs2[i] + dsum2_c6[i] * nz_mod) / (1.0 * ne[i])
    else:
        zeff_gph_c6[i] = dsum2_c6[i] / dsum1_c6[i]

# Need the densities for C6+ and all carbon to compare.
divimp_nz_all = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge="all")["nz"]) * gop_absfac * nz_mod
divimp_nz_c6 = np.array(gop.fake_probe(1.9, 2.5, 0.0, 0.0, "nz", charge=6)["nz"]) * gop_absfac * nz_mod

# Distill into simple arrays.
tchords_r = [cer[c]["r"] for c in cer.keys() if c[0]=="T"]
tchords_zeff = [cer[c]["zeff"] for c in cer.keys() if c[0]=="T"]
tchords_nz = [cer[c]["nz"] for c in cer.keys() if c[0]=="T"]
tchords_zeff_err = [cer[c]["zeff_err"] for c in cer.keys() if c[0]=="T"]
tchords_nz_err = [cer[c]["nz_err"] for c in cer.keys() if c[0]=="T"]
vchords_r = [cer[c]["r"] for c in cer.keys() if c[0]=="V"]
vchords_zeff = [cer[c]["zeff"] for c in cer.keys() if c[0]=="V"]
vchords_nz = [cer[c]["nz"] for c in cer.keys() if c[0]=="V"]
vchords_zeff_err = [cer[c]["zeff_err"] for c in cer.keys() if c[0]=="V"]
vchords_nz_err = [cer[c]["nz_err"] for c in cer.keys() if c[0]=="V"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4), sharex=True)

ax1.axvline(2.2455, color="k", linestyle="--")
ax1.scatter(tchords_r, tchords_zeff, c="tab:red", marker="o", label="T")
ax1.scatter(vchords_r, vchords_zeff, c="tab:red", marker="*", label="V")
ax1.plot(divimp_r, zeff_gph, color="tab:red", label="DIVIMP-Graphite")
ax1.plot(divimp_r, zeff_gph_c6, color="tab:red", label="DIVIMP-Graphite (C6+)", linestyle="--")
ax1.legend()
ax1.set_ylim([0.0, 3])
ax1.set_xlim([2.1, 2.35])
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Zeff")

ax2.axvline(2.2455, color="k", linestyle="--")
ax2.scatter(tchords_r, tchords_nz, c="tab:red", marker="o", label="T")
ax2.scatter(vchords_r, vchords_nz, c="tab:red", marker="*", label="V")
ax2.plot(divimp_r, divimp_nz_all, color="tab:red", label="DIVIMP-Graphite")
ax2.plot(divimp_r, divimp_nz_c6, color="tab:red", linestyle="--", label="DIVIMP-Graphite (C6+)")
ax2.legend()
ax2.set_xlabel("R (m)")
ax2.set_ylabel("Carbon Density (m-3)")

fig.suptitle("nz_mod = {}".format(nz_mod))
fig.tight_layout()
fig.show()
