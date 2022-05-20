# Script to choose from a select amount of shots and to compare the resulting
# dominant frequencies.
from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate
from matplotlib.colors import LogNorm


# Ideally LSN shots.
# L: 167196, 167200, 176508, 179900, 184948, 184952 (these two are odd), 164262, 167277
# H: 167247, 167277, 170848 (ramp)
shots = [167196, 167200, 170843, 170844, 170845, 170846, 176508, 179900, 164262, 167277]
#shots = [164245, 164253]

# Show the unsmoothed frequency spectrum. Best with just a single shot, specify
# below.
show_unsmooth = False
if show_unsmooth:
    shots = [167195]
    lowcut = 1500
    highcut = 2400

# Default: Restrict to a certain time range to avoid ramp up/down.
time_range = [2000, 4500]

# Tile to monitor.
#tilename = "ILEG15B000"
tilename = "ILEG15B150"

# First just load the data.
tiledata = {}
for shot in shots:
    print(shot)
    conn = MDSplus.Connection("atlas.gat.com")
    gaobj = gadata(tilename, shot, connection=conn, print_out=False)
    t = np.array(gaobj.xdata)
    y = np.array(gaobj.zdata)
    tiledata[shot] = {"t":t, "y":y}

# For each shot, perform a FFT to get the frequency data.
for shot in shots:
    t = tiledata[shot]["t"]
    y = tiledata[shot]["y"]

    # Per-shot time ranges if necessary.
    if shot in [170843, 170844, 170845, 170846]:
        mask = np.logical_and(t>2000, t<3800)
    elif shot in [164262, 164245, 164253]:
        #mask = np.logical_or(np.logical_and(t>1500, t<1900), np.logical_and(t>4700, t<5200))
        mask = np.logical_and(t>1500, t<1900)
    elif shot in [148673]:
        mask = np.logical_and(t>1500, t<2800)
    elif shot in [172409, 172413]:
        mask = np.logical_and(t>1600, t<2000)
    elif shot in [174212]:
        mask = np.logical_and(t>3000, t<4000)
    elif shot in [174305]:
        mask = np.logical_and(t>2500, t<3000)

    else:
        mask = np.logical_and(t>time_range[0], t<time_range[1])
    t = t[mask]
    y = y[mask]

    # Fourier transform. Calculate average sample rate (should be ~constant).
    #sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)]) * 1000  # 1/ms to 1/s
    sample_rate = 100000   # 100 kHz
    yf = rfft(y)
    xf = rfftfreq(len(t), 1 / sample_rate)

    tiledata[shot]["xf"] = xf
    tiledata[shot]["yf"] = np.abs(yf)
    tiledata[shot]["x"] = t
    tiledata[shot]["y"] = y

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot the frequencies on each other.
for shot in shots:
    if show_unsmooth:
        ax1.plot(tiledata[shot]["xf"], tiledata[shot]["yf"], label=shot, alpha=0.3)
    ax1.plot(tiledata[shot]["xf"], savgol_filter(tiledata[shot]["yf"], 501, 5), label=shot)
ax1.legend(loc="upper right")
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Magnitude")
ax1.set_yscale("log")

# Only do if show_unsmooth since that's typically just one shot which is fine here.
# Apply a bandpass filter.
shot = shots[0]
nyq = 100000 * 0.5
low = lowcut / nyq
high = highcut / nyq
#b, a = butter(2, high, btype='lowpass')
b, a = butter(2, [low, high], btype='bandpass')
y = filtfilt(b, a, tiledata[shot]["y"])

ax2.plot(tiledata[shot]["x"], y)
ax2.set_xlabel("Time (ms)")
ax2.set_ylabel("Amplitude")

fig.suptitle(tilename)
fig.tight_layout()
fig.show()
