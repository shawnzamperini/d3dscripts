from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate
from ThomsonClass import ThomsonClass


shot    = 176508
time_range = [2000, 5000]
lowcut  = 500
highcut = 5000


ts = ThomsonClass(shot, 'core')
ts.load_ts(tunnel=False)
ts_dict = ts.ts_dict
signals = range(0, 21)
for signal in signals:
    t = ts_dict["density"]["X"]
    y = ts_dict["density"]["Y"][signal]
    mask = np.logical_and(t>time_range[0], t<time_range[1])
    t = t[mask]
    y = y[mask]

    # Fourier transform. Calculate average sample rate (should be ~constant).
    sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)])
    print("  {:d} Hz".format(int(sample_rate)))
    yf = fft(y)
    xf = fftfreq(len(t), 1 / sample_rate)

    # Divide by the mean from 1000 Hz onwards as a way of normalizing.
    #avg = np.abs(yf[xf > 1000]).mean()

    # Ignore anything below zero Hz.
    mask = xf > 0
    xf = xf[mask]
    yf = yf[mask]

    # Apply a bandpass filter in the range of the blobs.
    nyq = sample_rate * 0.5
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(2, [low, high], btype='bandpass')
    yfilt = filtfilt(b, a, y)

    data[signal] = {"xf":xf, "yf":np.abs(yf), "rate":sample_rate,
        "x":t, "y":y, "yfilt":yfilt}


# Plotting.
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9,7))

for signal in good_signals:
    ax1.plot(data[signal]["xf"], savgol_filter(data[signal]["yf"], 201, 2), label=signal)
    ax2.plot(data[signal]["x"], data[signal]["yfilt"], label=signal)

#ax1.set_xlim([0, 5000])
#ax1.set_ylim([0, 10])
ax1.legend(loc="upper right")
ax2.legend(loc="upper right")

fig.tight_layout()
fig.show()
