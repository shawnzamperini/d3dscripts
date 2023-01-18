# This script take the radial filterscope signals and get that's radial derivative for use in the simple SOL equations
# as a radial flux of ionization.
import MDSplus
import numpy as np
import matplotlib.pyplot as plt
from gadata import gadata
import sys


# Unsure how to pull the R locations, so eyeballing them and manually entering them here.
chord_rs = [2.207, 2.220, 2.233, 2.246, 2.258, 2.269, 2.281, 2.293]

def main(shot, time, sxb=10):

    # Load, extract data for each chord.
    conn = MDSplus.Connection("atlas.gat.com")
    fs = {}
    for chord in range(1, 9):
        print("Loading chord {}...".format(chord))
        tag = "FS{}MIDDA".format(chord)
        gaobj = gadata(tag, shot, connection=conn)
        times = gaobj.xdata
        signal = 4 * np.pi * gaobj.zdata * sxb * 10000 # cm2 to m2

        # Get an average value around a 100 ms time window.
        mask = np.abs(times-time) > 100
        fs[chord] = signal[mask].mean()

    chord_data = list(fs.values())


    fig, ax1 = plt.subplots()
    ax1.scatter(chord_rs, chord_data)
    ax1.set_xlabel("R (m)")
    ax1.set_ylabel("Ionization flux (ions/m2/s)")
    fig.tight_layout()
    fig.show()


    return {"gaobj":gaobj, "fs":fs}

if __name__ == "__main__":
    shot = int(sys.argv[1])
    time = float(sys.argv[2])
    a = main(shot, time)
