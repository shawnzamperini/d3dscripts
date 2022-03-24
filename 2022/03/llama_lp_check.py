# This script looks at some of the LP data to see if maybe the peaks lines up
# with detahcment.
import get_lp
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from gadata import gadata
import MDSplus as mds


# Load LP data for a particular LP. Apply some smoothing.
#shot = 180910
#pname = "probe 127"
pname = "probe 19"

data = {}
for shot in [180910, 180914, 180916]:

    if shot == 180914:
        tmin = 2700
        tmax = 4400
    elif shot == 180916:
        tmin = 2700
        tmax = 4400
    elif shot == 180910:
        tmin = 2000
        tmax = 3500

    lp_dict = get_lp.get_dict_of_lps(shot, False)
    #lp_dict = get_lp.plot_lps(shot, tmin, tmax, tunnel=False, bins=25, xtype="time")
    #mask = np.array(lp_dict["pnames"]) == pname
    #lp_x = np.array(lp_dict["time"])[mask]
    #lp_y = np.array(lp_dict["jsat (A/cm2)"])[mask]

    lp_x = lp_dict[pname]["time"]
    lp_y = savgol_filter(lp_dict[pname]["jsat"], 91, 2)

    # Load line averaged density.
    conn = mds.Connection("atlas.gat.com")
    ga = gadata("DENSITY", shot, connection=conn)
    ne_x = ga.xdata
    ne_y = ga.zdata * 1e6  # cm3 to m3

    # Interpolate values onto common time axis.
    ts = np.arange(tmin, tmax, 5)

    f_ne = interp1d(ne_x, ne_y)
    f_jsat = interp1d(lp_x, lp_y)
    ne_int = f_ne(ts)
    jsat_int = f_jsat(ts)

    # Store in dictionary.
    data[shot] = {"ts":ts, "jsat":jsat_int, "ne":ne_int}

fig, ax = plt.subplots()
colors = ["tab:red", "tab:purple", "tab:cyan"]
count = 0
for shot, vals in data.items():
    ax.scatter(vals["ne"], vals["jsat"], s=20, c=colors[count], edgecolors="k", label=shot)
    count += 1
ax.set_xlabel("line-averaged ne (m-3)")
ax.set_ylabel("jsat (A/cm2)")
ax.legend()
ax.set_title("{} jsat rollover curve".format(pname))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
