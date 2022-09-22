# Estimate Er at the window-frame for that set of shots.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))

for p in [1, 2]:

    shot = 190485
    plunge = p
    window = 2.289


    label = "{}_{}".format(shot, plunge)
    path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP{}_{}.tab".format(shot, plunge)
    df = pd.read_csv(path, sep="\t")

    r = df["R(cm)"].values / 100

    # Method #1 by using floating potential values. # Linear fit the values +/-
    # around the window region to get Er at the window.
    te = df["Te(eV)"].values
    vp = (df["Vf1(V)"]+df["Vf2(V)"]).values / 2 + 3 * te

    give = 0.01
    mask = np.where(np.logical_and(r>window-give, r<window+give))
    z = np.polyfit(r[mask], vp[mask], 1)
    er1 = -z[0]
    print("Er1 = {:.2f}".format(er1))

    # Method #2 by assuming Er ~ -3 * dTe/dr.
    er2 = -3 * np.gradient(te, r)

    ax1.axvline(window, color="k", linestyle="--")
    ax1.scatter(r, vp)
    ax1.scatter(r[mask], vp[mask])
    ax1.set_xlabel("R (m)")
    ax1.set_ylabel("Vp (V)")

    ax2.scatter(r, er2, label="-3 * dTe/dr")
    ax2.plot(r[mask], np.full(len(mask[0]), er1), label="dVp/dr")
    ax2.legend()
    ax2.grid()
    ax2.set_ylim([-10, 750])

fig.suptitle(label)
fig.tight_layout()
fig.show()
