from gadata import gadata
import numpy as np
import MDSplus
from scipy.signal import find_peaks, savgol_filter
import matplotlib.pyplot as plt


shot = 190178
signal = "FS02UP"
#signal = "FS4MIDDA"
times = [2000, 5000]
width = 200
height = 0.3e16

conn = MDSplus.Connection("atlas.gat.com")
gaobj = gadata(signal, shot, connection=conn)
time = gaobj.xdata
fs = gaobj.zdata
mask = np.logical_and(time>times[0], time<times[1])
time = time[mask]
fs = fs[mask]

peaks = find_peaks(fs, distance=200, height=height)
peak_idx = peaks[0]

# Count number of peaks in each bin.
bins = np.arange(times[0], times[1], width)
counts = np.zeros(bins.shape)
for i in range(0, len(peak_idx)):
    t = time[peak_idx[i]]
    for j in range(0, len(bins)-1):
        if t >= bins[j] and t < bins[j+1]:
            counts[j] += 1
freqs = counts / (width/1000)
bin_mids = [bins[i] + (bins[i+1] - bins[i])/2 for i in range(0, len(bins)-1)]

# Load SPRED FE23 signal.
gaobj2 = gadata("SPRED_FE23", shot, connection=conn)
time2 = gaobj2.xdata
fe23 = savgol_filter(gaobj2.zdata, 31, 2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
ax1.plot(time, fs, color="k")
ax1.scatter(time[peak_idx], fs[peak_idx], s=10, color="r")
ax1.set_xlabel("Time (ms)")
ax1.set_ylabel(signal)
fig.tight_layout()
fig.show()

fig, ax2 = plt.subplots(figsize=(5,4))
ax2.plot(bin_mids, freqs[:-1], color="tab:red")
ax2.set_xlabel("Time (ms)", fontsize=14)
ax2.tick_params(axis="y", colors="tab:red")
ax2.set_ylabel("ELM Frequency (Hz)", color="tab:red", fontsize=14)
ax2.grid()
ax22 = ax2.twinx()
ax22.plot(time2, fe23, color="tab:purple")
ax22.set_ylabel("SPRED_Fe23", color="tab:purple", fontsize=14)
ax22.tick_params(axis='y', colors='tab:purple')
ax22.set_ylim([0, 6e13])
ax2.set_xlim(times)
fig.suptitle(shot)
fig.tight_layout()
fig.show()
