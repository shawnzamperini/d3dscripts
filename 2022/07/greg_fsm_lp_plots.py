# A script to plot some LP data against accompanying filterscope data.
import get_lp
import matplotlib.pyplot as plt
from gadata import gadata
import MDSplus
from scipy.interpolate import interp1d
import numpy as np
from scipy.signal import savgol_filter


shot = 190444  # Attached
tmin = 2500
tmax = 5500
filt_sig = "tube39:pmt_volt"  # VW2 W1M 39
lp_name = "probe 163"  # A12 for 190444

# Load lP data.
lp_dict = get_lp.plot_lps(shot, tmin, tmax, tunnel=False, xtype="time", bins=25)

# Load filterscope data.
gaobj = gadata(filt_sig, shot, tree="fscope", connection=MDSplus.Connection("atlas.gat.com"))

# Extract data.
mask = np.array(lp_dict["pnames"]) == lp_name
lp_t = np.array(lp_dict["time"])[mask]
lp_ne = np.array(lp_dict["ne (cm-3)"])[mask]
lp_te = np.array(lp_dict["Te (eV)"])[mask]
sort_idx = np.argsort(lp_t)
lp_t = lp_t[sort_idx]
lp_te = lp_te[sort_idx]
lp_ne = lp_ne[sort_idx]

filt_t = gaobj.xdata
filt_y = savgol_filter(gaobj.zdata, 1001, 2)

fig, ax1 = plt.subplots(figsize=(5,4))

#ax11 = ax1.twinx()

ax1.plot(filt_t, filt_y*500, color="k")
ax1.plot(lp_t, lp_te, color="r", lw=3)
ax1.plot(lp_t, lp_ne/1e12, color="g", lw=3)

ax1.set_xlim([2500, 5400])
ax1.set_ylim([0, None])
ax1.set_xlabel("Time (ms)", fontsize=16)
#ax1.set_ylabel("Filterscope W Signal", fontsize=16)
#ax11.set_ylabel("Langmuir Probe", fontsize=16)
ax1.set_title(shot, fontsize=16)
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)
ax1.tick_params(labelsize=12)

fig.tight_layout()
fig.show()
