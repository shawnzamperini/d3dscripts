import matplotlib.pyplot as plt
import numpy as np
from gadata import gadata
import MDSplus
import pickle

tchords = ["T{}".format(i) for i in range(1, 49)]
vchords = ["V{}".format(i) for i in range(1, 33)]
chords = np.append(tchords, vchords)

conn = MDSplus.Connection("atlas.gat.com")
shot = 163150
tmin = 1400
tmax = 1600
cer = {}
for chord in chords:
    print("Chord: {}".format(chord))
    signal_r = "CERAR{}".format(chord)
    signal_z = "CERAZ{}".format(chord)
    signal_zeff = "CERAZEFF{}".format(chord)
    signal_nz = "CERANZ{}".format(chord)
    signal_tz = "CERATI{}".format(chord)
    gaobj_r = gadata(signal_r, shot, connection=conn)
    gaobj_z = gadata(signal_z, shot, connection=conn)
    gaobj_zeff = gadata(signal_zeff, shot, connection=conn)
    gaobj_nz = gadata(signal_nz, shot, connection=conn)
    gaobj_tz = gadata(signal_tz, shot, connection=conn)
    chord_r = float(gaobj_r.zdata[0])
    chord_z = float(gaobj_z.zdata[0])

    # Return average Zeff value.
    mask = np.logical_and(gaobj_zeff.xdata > tmin, gaobj_zeff.xdata < tmax)
    chord_zeff = float(gaobj_zeff.zdata[mask].mean())
    chord_nz = float(gaobj_nz.zdata[mask].mean())
    chord_zeff_std = float(gaobj_zeff.zdata[mask].std())
    chord_nz_std = float(gaobj_nz.zdata[mask].std())

    # The temperature has its own time series.
    mask = np.logical_and(gaobj_tz.xdata > tmin, gaobj_tz.xdata < tmax)
    chord_tz = float(gaobj_tz.zdata[mask].mean())
    chord_tz_std = float(gaobj_tz.zdata[mask].std())

    # If zero, then no data for this chord.
    if chord_zeff != 0 and ~np.isnan(chord_zeff):
        cer[chord] = {"r": chord_r, "z": chord_z, "zeff": chord_zeff, "tz": chord_tz, "tz_err": chord_tz_std,
                      "nz": chord_nz, "zeff_err": chord_zeff_std, "nz_err": chord_nz_std}

print("Data found for following chords: {}".format(cer.keys()))

# Pickle data and save.
fname = "cer_{}_{}_{}.pickle".format(shot, tmin, tmax)
with open(fname, "wb") as f:
    pickle.dump(cer, f)

