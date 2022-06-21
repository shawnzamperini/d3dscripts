# Do FFT's of a broad overview of the MRC shot.
from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, butter, filtfilt


# Choose the type of shots with the most data: H-mode, shelf SP.
shots = [167247, 167276, 167277, 167279, 167320, 167528]
times = {167247:[2500, 5000],
167276:[2500, 5000],
167277:[2500, 4000],
167279:[2500, 5000],
167320:[2500, 5000]}
conn = MDSplus.Connection("atlas.gat.com")
signal = "ILEG15B150"

# Go through one shot at a time, load the tile data, perform a FFT on it and
# then save.
data = {}
for shot in shots:
    print(shot)
    gaobj = gadata(signal, shot, connection=conn)
    t = np.array(gaobj.xdata) / 1000
    y = np.array(gaobj.zdata)
    if shot == 167528:
        mask = np.logical_and(t>1.5, t<2.9)
        t1 = t[mask]
        y1 = y[mask]
        sample_rate = 100000   # 100 kHz
        yf = rfft(y1)
        xf = rfftfreq(len(t1), 1 / sample_rate)
        data["167528_1"] = {}
        data["167528_1"]["xf"] = xf
        data["167528_1"]["yf"] = np.abs(yf)
        data["167528_1"]["x"] = t1
        data["167528_1"]["y"] = y1

        mask = np.logical_and(t>3.5, t<5.0)
        t2 = t[mask]
        y2 = y[mask]
        sample_rate = 100000   # 100 kHz
        yf = rfft(y2)
        xf = rfftfreq(len(t2), 1 / sample_rate)
        data["167528_2"] = {}
        data["167528_2"]["xf"] = xf
        data["167528_2"]["yf"] = np.abs(yf)
        data["167528_2"]["x"] = t2
        data["167528_2"]["y"] = y2
    else:
        mask = np.logical_and(t>times[shot][0]/1000, t<times[shot][1]/1000)
        t = t[mask]
        y = y[mask]

        sample_rate = 100000   # 100 kHz
        yf = rfft(y)
        xf = rfftfreq(len(t), 1 / sample_rate)

        data[shot] = {}
        data[shot]["xf"] = xf
        data[shot]["yf"] = np.abs(yf)
        data[shot]["x"] = t
        data[shot]["y"] = y

# Plot 5 shots per plot.
fig, ax1 = plt.subplots(sharex=True, sharey=True)

for shot in shots[:5]:
    window = 201
    yfilt = savgol_filter(data[shot]["yf"], window, 2)
    ax1.plot(data[shot]["xf"], yfilt, label=shot)

ax1.legend()
ax1.set_xlim([0, 6000])
ax1.set_ylim([0, 150000])
fig.tight_layout()
fig.show()

nyq = 100000 * 0.5
#low = lowcut / nyq
high = 7500 / nyq
b, a = butter(2, high, btype='lowpass')
#b, a = butter(2, [low, high], btype='bandpass')
y = filtfilt(b, a, data["167528_1"]["y"])

fig, (ax1, ax2) = plt.subplots(1, 2)
window = 201
yfilt1 = savgol_filter(data["167528_1"]["yf"], window, 2)
yfilt2 = savgol_filter(data["167528_2"]["yf"], window, 2)
ax1.plot(data["167528_1"]["xf"], yfilt1, label="Ohmic")
ax1.plot(data["167528_2"]["xf"], yfilt2, label="H-mode")
ax1.legend()
ax2.plot(data["167528_1"]["x"], y)
fig.tight_layout()
fig.show()
