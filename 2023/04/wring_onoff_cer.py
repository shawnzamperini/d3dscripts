import MDSplus
import numpy as np
import matplotlib.pyplot as plt
from gadata import gadata



tchords = ["T{}".format(i) for i in range(1, 49)]
vchords = ["V{}".format(i) for i in range(1, 33)]
chords = np.append(tchords, vchords)

conn = MDSplus.Connection("atlas.gat.com")
tmin = 2500
tmax = 5000
shot_data = {}
for shot in [167190, 167193]:
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
    shot_data[shot] = cer

# Distill into simple arrays and plot.
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
fig, ax1 = plt.subplots(figsize=(5, 4))
for shot in shot_data.keys():
    cer = shot_data[shot]
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

    if shot == 167190:
        color = "k"
        label = "Graphite"
    elif shot == 167193:
        color = "hotpink"
        label = "Tungsten"
    ax1.errorbar(tchords_r, tchords_nz, yerr=tchords_nz_err, ms=5, color=color, markeredgecolor="k", lw=0, marker="o",
                 elinewidth=2, label=label)
    # ax2.errorbar(tchords_r, tchords_zeff, yerr=tchords_zeff_err, ms=5, color=color, markeredgecolor="k", lw=0, marker="o",
    #              elinewidth=2)

ax1.set_xlim(1.6, None)
ax1.set_xlabel("R (cm)", fontsize=14)
ax1.set_ylabel(r"Carbon Density $(\mathdefault{m^{-3}})$", fontsize=14)
ax1.legend(fontsize=14)
fig.tight_layout()
fig.show()
