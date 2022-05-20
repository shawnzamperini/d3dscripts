# Do a FFT on a filterscope signal.
from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter
from matplotlib.colors import LogNorm


# Inputs
shot = 167196
signal = "FS06"
time_range = [2, 5]
measured_fblob = 1860  # Hz
cbar_lims = [1e20, 1e23]

if signal == "FS06":
    ylims = [0, 3e18]
    vs = [1e23, 1e26]

# Grab the data, restrict to range of interest.
conn = MDSplus.Connection("atlas.gat.com")
gaobj = gadata(signal, shot, connection=conn)
t = gaobj.xdata / 1000
y = gaobj.zdata
mask = np.logical_and(t>time_range[0], t<time_range[1])
t = t[mask]
y = y[mask]

# Fourier transform. Calculate average sample rate (should be ~constant).
sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)])
yf = fft(y)
xf = fftfreq(len(t), 1 / sample_rate)

# Create spectrogram for plotting.
fs, ts, Sxx = spectrogram(y, fs=sample_rate)

# Plotting.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4))
ax1.plot(xf, np.abs(yf), color="k")
ax1.set_xlim([0, 1e4])
ax1.set_ylim(ylims)
#ax1.set_yscale("log")
ax1.set_xlabel("Frequency (Hz)")
ax1.axvline(measured_fblob, color="k", linestyle="--")
ax1.plot(xf, savgol_filter(np.abs(yf), 101, 2), color="r")

# A spectrogram.
cmesh = ax2.pcolormesh(ts+t.min(), fs, Sxx, shading="gouraud", norm=LogNorm(vmin=vs[0], vmax=vs[1]), cmap="inferno")
#cmesh = ax2.pcolormesh(ts+t.min(), fs, Sxx, shading="gouraud", norm=LogNorm(), cmap="inferno")
cbar = fig.colorbar(cmesh, ax=ax2)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Frequency (Hz)")
cbar.set_label("Magnitude (?)")
ax2.axhline(measured_fblob, color="k", linestyle="--")
#ax2.set_yscale("log")
ax2.set_ylim([0, 1e4])

fig.tight_layout()
fig.show()
