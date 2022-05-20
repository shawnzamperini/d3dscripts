# Do a Forier transform on the raw Langmuir probe data.
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt, savgol_filter
from scipy.fft import rfft, rfftfreq
import MDSplus
from gadata import gadata


# Inputs.
shot = 167196
signal = "TPLANG50"
low = 800
high = 1200
plot_range = [2498, 2502]


# Load data.
conn = MDSplus.Connection("atlas.gat.com")
gaobj = gadata(signal, shot, connection=conn)
time = np.array(gaobj.xdata)
curr = np.array(gaobj.zdata)

# Perform Fourier transform.
yf = np.abs(rfft(curr))
xf = rfftfreq(len(time), 1/1e6)

# Perform a bandstop filter.
nyq = 1e6 * 0.5
b, a = butter(2, [low/nyq, high/nyq], btype="bandstop")
yfilt = filtfilt(b, a, curr)

# Plotting.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(xf, (yf-yf.mean())/yf.std())
ax1.axvline(low, color="k")
ax1.axvline(high, color="k")
ax1.set_xlim([0, 5000])
ax1.set_ylim([0, 2])
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Magnitude (Arbitrary)")

ax2.plot(time, (curr-curr.mean())/curr.std(), color="k")
ax2.plot(time, (yfilt-yfilt.mean())/yfilt.std(), color="r")
ax2.set_xlabel("Time (ms)")
ax2.set_ylabel("Current (A?)")
ax2.set_xlim(plot_range)

fig.tight_layout()
fig.show()
