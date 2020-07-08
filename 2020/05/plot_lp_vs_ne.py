import get_lp
from gadata import gadata
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from scipy.signal import savgol_filter


# Grab the LP data and pedestal density.
d = get_lp.get_dict_of_lps(167247)
neped = gadata('prmtan_neped', 167247)

def interp_func(x, y):

    # Strip zeros.
    keep_idx = x != 0
    x = x[keep_idx]
    y = y[keep_idx]
    return interp1d(x, y)

# Interpolation functions so we can use the same time ranges.
f11 = interp_func(d['probe 11']['time'], d['probe 11']['jsat'])
f13 = interp_func(d['probe 13']['time'], d['probe 13']['jsat'])
f15 = interp_func(d['probe 15']['time'], d['probe 15']['jsat'])
f17 = interp_func(d['probe 17']['time'], d['probe 17']['jsat'])
f19 = interp_func(d['probe 19']['time'], d['probe 19']['jsat'])
f21 = interp_func(d['probe 21']['time'], d['probe 21']['jsat'])
fneped = interp_func(neped.xdata, neped.zdata)

# A common time range.
tcom = np.linspace(1400, 1600, 125)

# Sort according to neped (okay only if neped constantly increasing).
sort_idx = np.argsort(fneped(tcom))

fig, ax = plt.subplots()
#ax.plot(savgol_filter(fneped(tcom)[sort_idx], 11, 3), savgol_filter(f11(tcom)[sort_idx], 11, 3), label='11')
ax.plot(savgol_filter(fneped(tcom)[sort_idx], 17, 3), savgol_filter(f13(tcom)[sort_idx], 17, 3), label='13')
ax.plot(savgol_filter(fneped(tcom)[sort_idx], 17, 3), savgol_filter(f15(tcom)[sort_idx], 17, 3), label='15')
ax.plot(savgol_filter(fneped(tcom)[sort_idx], 17, 3), savgol_filter(f17(tcom)[sort_idx], 17, 3), label='17')
ax.plot(savgol_filter(fneped(tcom)[sort_idx], 17, 3), savgol_filter(f19(tcom)[sort_idx], 17, 3), label='19')
ax.plot(savgol_filter(fneped(tcom)[sort_idx], 17, 3), savgol_filter(f21(tcom)[sort_idx], 17, 3), label='21', lw=5)
ax.set_xlabel('neped (m-3)', fontsize=16)
ax.set_ylabel('jsat (A/cm2)', fontsize=16)
ax.legend(fontsize=16)
fig.tight_layout()
fig.show()
