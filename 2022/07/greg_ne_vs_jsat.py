import get_lp
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from gadata import gadata
import MDSplus



# Load data up front.
shots = [190485]
colors=["k", "r"]
tstart = 2800
tend = 5500
probe = "probe 161"   # MDS index 73 7/7/22
fig, ax1 = plt.subplots(figsize=(5, 4))
for i in range(0, len(shots)):

    shot = shots[i]
    lp_dict = get_lp.get_dict_of_lps(shot, tunnel=False)
    gaobj = gadata("DENSV2", shot, connection=MDSplus.Connection("atlas.gat.com"))

    # Pull out particular LP.
    lpt = lp_dict[probe]

    # Create interpolation objects.
    fne = interp1d(gaobj.xdata, gaobj.zdata)
    fjsat = interp1d(lp_dict[probe]["time"], lp_dict[probe]["jsat"])

    tmin = max(gaobj.xdata.min(), min(lp_dict[probe]["time"]), tstart)
    tmax = min(gaobj.xdata.max(), max(lp_dict[probe]["time"]), tend)

    time = np.linspace(tmin, tmax, 400)
    ne = fne(time)
    jsat = fjsat(time)


    ax1.scatter(ne, jsat, color=colors[i], alpha=0.75, marker="^", s=100, zorder=10)
ax1.grid(zorder=1)
ax1.set_title(shot, fontsize=20)
ax1.set_xlabel("ne (m-3)", fontsize=16)
ax1.set_ylabel("jsat (A/cm2)", fontsize=16)
ax1.set_ylim([0, 30])
fig.tight_layout()
fig.show()
