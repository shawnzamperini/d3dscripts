from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate
from matplotlib.colors import LogNorm


shot = 190060
signal = "FS4MIDDA"
time_range = [2, 3.5]

data = {}
conn = MDSplus.Connection("atlas.gat.com")
gaobj = gadata(signal, shot, connection=conn)
t = np.array(gaobj.xdata) / 1000
y = np.array(gaobj.zdata)
mask = np.logical_and(t>time_range[0], t<time_range[1])
t = t[mask]
y = y[mask]

# Fourier transform. Calculate average sample rate (should be ~constant).
#sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)])
sample_rate = 100000
yf = rfft(y)
xf = rfftfreq(len(t), 1 / sample_rate)

# Divide by the mean from 1000 Hz onwards as a way of normalizing.
avg = np.abs(yf[xf > 1000]).mean()

# Ignore anything below zero Hz.
mask = xf > 0
xf = xf[mask]
yf = yf[mask]
yf = np.abs(yf)

data[signal] = {"xf":xf, "yf":np.abs(yf) / avg, "rate":sample_rate, "x":t, "y":y}

fig, ax1 = plt.subplots()

ax1.plot(xf, yf, alpha=0.2, color="k")
ax1.plot(xf, savgol_filter(yf, 31, 2), color="k")

time_range = [3.9, 5.0]
gaobj = gadata(signal, shot, connection=conn)
t = np.array(gaobj.xdata) / 1000
y = np.array(gaobj.zdata)
mask = np.logical_and(t>time_range[0], t<time_range[1])
t = t[mask]
y = y[mask]

# Fourier transform. Calculate average sample rate (should be ~constant).
#sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)])
sample_rate = 100000
yf = rfft(y)
xf = rfftfreq(len(t), 1 / sample_rate)

# Divide by the mean from 1000 Hz onwards as a way of normalizing.
avg = np.abs(yf[xf > 1000]).mean()

# Ignore anything below zero Hz.
mask = xf > 0
xf = xf[mask]
yf = yf[mask]
yf = np.abs(yf)

ax1.plot(xf, yf, alpha=0.2, color="r")
ax1.plot(xf, savgol_filter(yf, 31, 2), color="r")

ax1.set_xlabel("Frequency")
ax1.set_ylabel("Amplitude")

fig.tight_layout()
fig.show()
