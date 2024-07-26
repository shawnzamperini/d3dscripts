# Is there any trend in the Fe23 signal from SPRED with q95 (indicative of
# tile misalignments)?
from gadata import gadata
import numpy as np
import matplotlib.pyplot as plt
import MDSplus
from scipy.interpolate import interp1d
from scipy.signal import medfilt


shots = [190448, 190449, 190450, 190451, 190452]
tmin = 1580
tmax = 4000

results = {}
conn = MDSplus.Connection("atlas.gat.com")
for shot in shots:

    # Load signals.
    print(shot)
    gaobj = gadata("SPRED_FE23", shot, connection=conn)
    fe23_t = gaobj.xdata
    fe23 = gaobj.zdata
    gaobj = gadata("Q95", shot, connection=conn)
    q95_t = gaobj.xdata
    q95 = gaobj.zdata

    # Need to smooth fe23.
    fe23s = medfilt(fe23, 21)

    # Need to interpolate onto a common time axis.
    f_q95 = interp1d(q95_t, q95)

    # Get values of each just for the time range.
    fe23_mask = np.logical_and(fe23_t>tmin, fe23_t<tmax)
    fe23_keep = fe23[fe23_mask]
    fe23s_keep = fe23s[fe23_mask]
    q95_keep = f_q95(fe23_t[fe23_mask])

    # Find max fe23 and the corresponding q95 value.
    max_idx = np.argmax(fe23s_keep)
    fe23_max = fe23s_keep[max_idx]
    q95_max = q95_keep[max_idx]

    results[shot] = {"q95":q95_keep, "fe23":fe23_keep, "fe23s":fe23s_keep,
        "fe23_max":fe23_max, "q95_max":q95_max, "t":fe23_t[fe23_mask],
        "max_idx":max_idx}

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))

colors = ["C{:}".format(i) for i in range(0, len(shots))]
count = 0
for shot, data in results.items():
    ax1.plot(data["t"], data["fe23"], alpha=0.3, color=colors[count], zorder=5)
    ax1.plot(data["t"], data["fe23s"], color=colors[count], zorder=10)
    ax1.scatter(data["t"][data["max_idx"]], data["fe23s"][data["max_idx"]], color=colors[count], marker="*", edgecolor="k", zorder=15)
    ax2.plot(data["t"], data["q95"], color=colors[count])
    ax3.scatter(data["q95_max"], data["fe23_max"], label=shot, color=colors[count])
    count += 1

ax3.legend()
ax3.set_xlabel("q95")
ax3.set_ylabel("fe23")
fig.tight_layout()
fig.show()
