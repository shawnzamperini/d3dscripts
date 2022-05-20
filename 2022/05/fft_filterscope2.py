from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate
from matplotlib.colors import LogNorm


# Inputs
shot = 172410
signals = ["FS0{}".format(i) for i in range(0, 9)]
#signals = np.append(signals, ["FS{}MIDDA".format(i) for i in range(1,9)])
#signals = np.append(signals, ["FS01DWDA", "FS04DWDA", "FS06DWDA", "FS11DWDA", "FS09DWDA", "FS12"])
signals = ["ILEG15B000","ILEG13B008","ILEG15B060","ILEG13B068","ILEG15B150","ILEG13B158","ILEG15B240","ILEG13B248","ILEG15B315","ILEG13B315"]
time_range = [2, 4]
measured_fblob = 1600  # Hz
cbar_lims = [1e20, 1e23]

# Grab the data, restrict to range of interest.
data = {}; good_signals = []
for signal in signals:
    print(signal)
    try:
        conn = MDSplus.Connection("atlas.gat.com")
        gaobj = gadata(signal, shot, connection=conn)
        t = np.array(gaobj.xdata) / 1000
        y = np.array(gaobj.zdata)
        mask = np.logical_and(t>time_range[0], t<time_range[1])
        t = t[mask]
        y = y[mask]

        # Fourier transform. Calculate average sample rate (should be ~constant).
        sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)])
        yf = fft(y)
        xf = fftfreq(len(t), 1 / sample_rate)

        # Divide by the mean from 1000 Hz onwards as a way of normalizing.
        avg = np.abs(yf[xf > 1000]).mean()

        # Ignore anything below zero Hz.
        mask = xf > 0
        xf = xf[mask]
        yf = yf[mask]

        data[signal] = {"xf":xf, "yf":np.abs(yf) / avg, "rate":sample_rate, "x":t, "y":y}
        good_signals.append(signal)
    except:
        print("  Skipped!")

# Plot them all smoothed.
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8,8), sharex=True, sharey=True)
for signal in good_signals:
    if signal[-5:] == "MIDDA":
        ax2.plot(data[signal]["xf"], savgol_filter(data[signal]["yf"], 201, 2), label=signal)
    elif signal in ["FS01DWDA", "FS04DWDA", "FS06DWDA", "FS11DWDA", "FS09DWDA", "FS12"]:
        ax3.plot(data[signal]["xf"], savgol_filter(data[signal]["yf"], 201, 2), label=signal)
    else:
        ax1.plot(data[signal]["xf"], savgol_filter(data[signal]["yf"], 501, 2), label=signal)
ax1.legend()
ax2.legend()
ax3.legend()
#ax1.set_xlim([0, 5000])
#ax1.set_ylim([0, 10])
fig.supylabel("Magnitude (arbitrary)")
fig.supxlabel("Frequency (Hz)")
#ax1.set_ylabel("Magnitude (arbitrary)")
#ax1.set_xlabel("Frequency")
ax1.axvline(measured_fblob, color="k", linestyle="--")
fig.tight_layout()
fig.show()

# Signal processing of some select signals. Copying example here:
# https://stackoverflow.com/questions/58846626/filter-out-range-of-frequencies-using-band-stop-filter-in-python-and-confirm-it
fig, ax1 = plt.subplots()
#for signal in ["FS03", "FS06", "FS1MIDDA", "FS4MIDDA"]:
#for signal in ["FS03"]:
#for signal in ["FS1MIDDA", "FS2MIDDA", "FS3MIDDA", "FS4MIDDA"]:
for signal in good_signals:
    nyq = data[signal]["rate"] * 0.5
    if signal in ["FS03", "FS06"]:
        lowcut = 1000
        highcut = 2000
    elif signal in ["FS1MIDDA", "FS4MIDDA"]:
        lowcut = 2000
        highcut = 3000
    else:
        lowcut = 2100
        highcut = 2500
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(2, high, btype='lowpass')
    #b, a = butter(2, [low, high], btype='bandpass')
    y = filtfilt(b, a, data[signal]["y"])

    #ax1.plot(data[signal]["x"], data[signal]["y"], alpha=0.3, color="r")
    ax1.plot(data[signal]["x"], y, label=signal)
    #ax1.plot(data[signal]["x"], savgol_filter(y, 201, 3))
    data[signal]["y_filt"] = y
ax1.legend(loc="upper right")
fig.tight_layout()
fig.show()

# Copying this function from https://stackoverflow.com/questions/41492882/find-time-shift-of-two-signals-using-cross-correlation
def lag_finder(y1, y2, sr):
    n = len(y1)

    corr = correlate(y2, y1, mode='same') / np.sqrt(correlate(y1, y1, mode='same')[int(n/2)] * correlate(y2, y2, mode='same')[int(n/2)])

    delay_arr = np.linspace(-0.5*n/sr, 0.5*n/sr, n)
    delay = delay_arr[np.argmax(corr)]
    print('y2 is ' + str(delay) + ' behind y1')

    plt.figure()
    plt.plot(delay_arr, corr)
    plt.title('Lag: ' + str(np.round(delay, 3)) + ' s')
    plt.xlabel('Lag')
    plt.ylabel('Correlation coeff')
    plt.show()

    return delay

delay = lag_finder(data["FS2MIDDA"]["y_filt"], data["FS4MIDDA"]["y_filt"], 50000)

# Then blob velocity is distance between the two location divided by the delay.
dist_btw = 2.246 - 2.221  # Eyeball approximation
vblob = dist_btw / delay
print("Estimated blob velocity = {:.2f}".format(vblob))

fig, ax1 = plt.subplots()
ax1.scatter(data["FS2MIDDA"]["y_filt"], data["FS4MIDDA"]["y_filt"])
fig.tight_layout()
fig.show()
