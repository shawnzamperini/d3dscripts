import matplotlib.pyplot as plt
from gadata import gadata
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


shot = 167350

def get_y_data(tag):
    ga_obj = gadata(tag, shot)
    return ga_obj.xdata, ga_obj.zdata

t_pinj, pinj = get_y_data("PINJ")
t_prad, prad = get_y_data("PRAD_CORE")

# Convert to MW.
pinj = pinj / 1000
prad = prad / 1000000

# Smooth PINJ.
pinj = savgol_filter(pinj, 31, 1)

# Interpolate onto the same domain.
f_pinj = interp1d(t_pinj, pinj)
pinj_int = f_pinj(t_prad)

psol = pinj_int - prad

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7), sharex=True, sharey=True)
ax1.plot(t_pinj, pinj, color="tab:red", lw=3)
ax2.plot(t_prad, prad, color="tab:red", lw=3)
ax1.set_ylabel("PINJ (MW)", fontsize=16)
ax2.set_ylabel("PRAD (MW)", fontsize=16)
ax2.set_xlabel("Time (ms)", fontsize=16)
ax1.set_xlim([0, 6000])
ax1.set_ylim([0, 3.5])
ax1.set_yticks([0, 1, 2, 3])
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax1.tick_params(which='both', labelsize=12)
ax2.tick_params(which='both', labelsize=12)
fig.tight_layout()
fig.show()

fig, ax = plt.subplots()
ax.plot(t_prad, psol, color="tab:red", lw=3)
ax.set_xlim([0, 6000])
ax.set_ylim([0, 3.5])
ax.set_yticks([0, 1, 2, 3])
ax.set_ylabel("PSOL (MW)", fontsize=16)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(which='both', labelsize=12)
fig.tight_layout()
fig.show()
