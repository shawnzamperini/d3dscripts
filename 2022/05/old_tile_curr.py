from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate
from matplotlib.colors import LogNorm


# Compare a shot that went form L to H mode and see how the tile current change.
shot = 119007
signal = "HCD34"
ltime = [1200, 2600]
htime = [2600, 3500]
stime = [4000, 4900]  # s for "sweep" time, which was after the strike point moved right up to the shelf, and went back into L-mode.
bandpass = [1500, 3000]

# Load data.
conn = MDSplus.Connection("atlas.gat.com")
gaobj = gadata(signal, shot, connection=conn)
time = np.array(gaobj.xdata)
tile = np.array(gaobj.zdata)
lmode = np.logical_and(time>ltime[0], time<ltime[1])
hmode = np.logical_and(time>htime[0], time<htime[1])
smode = np.logical_and(time>stime[0], time<stime[1])

# FFT.
sample_rate = 200000   # 200 kHz.
lyf = np.abs(rfft(tile[lmode]))
lxf = rfftfreq(len(time[lmode]), 1 / sample_rate)
hyf = np.abs(rfft(tile[hmode]))
hxf = rfftfreq(len(time[hmode]), 1 / sample_rate)
syf = np.abs(rfft(tile[smode]))
sxf = rfftfreq(len(time[smode]), 1 / sample_rate)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(lxf, savgol_filter(lyf, 201, 2), label="L-mode (floor)")
ax1.plot(hxf, savgol_filter(hyf, 201, 2), label="H-mode (floor)")
#ax1.plot(sxf, savgol_filter(syf, 501, 2), label="L-mode (shelf)")
ax1.plot(sxf, savgol_filter(syf, 21, 2), label="L-mode (shelf)")
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Spectral Intensity")
ax1.set_xlim([0, 6000])
ax1.set_ylim([0, 10000])
ax1.legend()

# Bandpass filtering of the range of interest.
nyq = sample_rate * 0.5
low = bandpass[0] / nyq
high = bandpass[1] / nyq
#b, a = butter(2, high, btype='lowpass')
b, a = butter(2, [low, high], btype='bandpass')
tile_filt = filtfilt(b, a, tile)
ax2.plot(time[lmode], tile_filt[lmode])
ax2.plot(time[hmode], tile_filt[hmode])
ax2.plot(time[smode], tile_filt[smode])
ax2.set_xlim([1200, 4800])
ax2.set_ylim([-7, 7])
ax2.set_ylabel("Tile Current {} Hz".format(bandpass))
ax2.set_xlabel("Time (ms)")
fig.tight_layout()
fig.show()

# Plot a representative inter-ELM section of the tile currents on some L-mode.
lwindow = [1840, 1880]
hwindow = [2840, 2880]
lsamp = tile_filt[np.logical_and(time>lwindow[0], time<lwindow[1])]
hsamp = tile_filt[np.logical_and(time>hwindow[0], time<hwindow[1])]
ltsamp = np.arange(0, len(lsamp)) / sample_rate
htsamp = np.arange(0, len(hsamp)) / sample_rate
fig, ax1 = plt.subplots()
ax1.plot(ltsamp, lsamp, label="L-mode")
ax1.plot(htsamp, hsamp, label="H-mode")
ax1.legend()
ax1.set_xlabel("Time (ms)")
ax1.set_ylabel("Tile Current {} Hz".format(bandpass))
fig.tight_layout()
fig.show()
