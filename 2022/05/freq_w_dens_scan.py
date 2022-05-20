# Purpose here is to look at how the peak frequency changes during a density scan.
from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate
from matplotlib.colors import LogNorm


#shots = [172410]
shots = [167195]

# Default: Restrict to a certain time range to avoid ramp up/down.
#time_range = [1600, 4000]
time_range = [3900, 5000]

# Tile to monitor.
#tilename = "ILEG15B000"
#tilename = "ILEG15B150"
tilename = "ILEG13B008"

# First just load the data.
tiledata = {}
for shot in shots:
    print(shot)
    conn = MDSplus.Connection("atlas.gat.com")
    gaobj = gadata(tilename, shot, connection=conn, print_out=False)
    t = np.array(gaobj.xdata)
    y = np.array(gaobj.zdata)
    gaobj = gadata("RVSOUT", shot, connection=conn, print_out=False)
    rvsout_t = np.array(gaobj.xdata)
    rvsout = np.array(gaobj.zdata)
    tiledata[shot] = {"t":t, "y":y, "rvsout_t":rvsout_t, "rvsout":rvsout}

# Time bins to average results in. Minimum bin value.
bins = np.linspace(time_range[0], time_range[1], 15)
avg_rvsouts = []

# Store a FFT for each time bin.
for shot in shots:
    t = tiledata[shot]["t"]
    y = tiledata[shot]["y"]
    rvsout_t = tiledata[shot]["rvsout_t"]
    rvsout = tiledata[shot]["rvsout"]

    for i in range(0, len(bins)):
        if i == len(bins)-1:
            continue
        else:
            mask = np.logical_and(t>=bins[i], t<bins[i+1])

        # Fourier transform.
        sample_rate = 100000   # 100 kHz
        yf = rfft(y[mask])
        xf = rfftfreq(len(t[mask]), 1 / sample_rate)

        tiledata[shot]["xf{}".format(i)] = xf
        tiledata[shot]["yf{}".format(i)] = np.abs(yf)

        if i == len(bins)-1:
            continue
        else:
            mask = np.logical_and(rvsout_t>=bins[i], rvsout_t<bins[i+1])
        avg_rvsouts.append(rvsout[mask].mean())


norm = plt.Normalize(bins.min(), bins.max())
cmap = plt.get_cmap("inferno")
fig, ax1 = plt.subplots(figsize=(11,5))
print("Max Frequency")
for i in range(0, len(bins)-1):
    x = tiledata[shots[0]]["xf{}".format(i)]
    y = tiledata[shots[0]]["yf{}".format(i)]
    ax1.plot(x, savgol_filter(y, 201, 2), color=cmap(norm(bins[i])), label=int(bins[i]))
    #ax1.plot(x, y, color=cmap(norm(bins[i])))

    # Ad-hoc way of pulling out peak frequency for 167195.
    window = [1280, 4000]
    mask = np.logical_and(x>window[0], x<window[1])
    idx = np.argmin(savgol_filter(y, 201, 2)[mask])
    print("{}".format(int(x[idx])))
    #print((bins[i] + bins[i+1])/2)

ax1.legend()
ax1.set_xlabel("Frequency")
ax1.set_ylabel("Magnitude (arbitrary)")
fig.tight_layout()
fig.show()

# Just plot a single time.
fig, ax1 = plt.subplots()
for i in [0, 7, 13]:
    x = tiledata[shots[0]]["xf{}".format(i)]
    y = tiledata[shots[0]]["yf{}".format(i)]
    #ax1.plot(x, y, color=cmap(norm(bins[i])), label=int(bins[i]), alpha=0.3)
    ax1.plot(x, savgol_filter(y, 201, 2), color=cmap(norm(bins[i])), label=int(bins[i]), lw=2)
ax1.set_xlim([0, 10000])
ax1.set_ylim([0, 10000])
#ax1.legend()
ax1.set_xlabel("Frequency", fontsize=16)
ax1.set_ylabel("Magnitude (arbitrary)", fontsize=16)
fig.tight_layout()
fig.show()
