# This script uses the midplane filterscope system to measure the radial
# blob velocity.
from gadata import gadata
import MDSplus
import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, savgol_filter, butter, filtfilt, correlate


# Inputs
shot = 186968
mode = "lower"  # One of lower, upper or midplane.

# Saved settings.
if shot == 184527:     # 50 kHz, methane experiment
    time_range = [2, 5]
    lowcut     = 300
    highcut    = 2000
    sep_chord  = "FS5MIDDA"     # Which filterscope chord to use for the "starting location" for the blobs.
elif shot == 176508:   # 100 kHz, SiC DiMES experiment
    time_range = [2, 5]
    lowcut     = 300
    highcut    = 2000
    sep_chord  = "FS2MIDDA"
    print("Warning: No real blob component in signal!")
elif shot in [167195, 167200]:   # 50 kHz, MRC classic
    time_range = [2, 5]
    lowcut     = 1000
    highcut    = 2000
    sep_chord  = "FS2MIDDA"
elif shot == 179900:   # 100 kHz, Si/SiC powder injection experiment. No significant blob component either...
    time_range = [2, 5]
    lowcut     = 1000
    highcut    = 2000
    sep_chord  = "FS2MIDDA"
    print("Warning: No real blob component in signal!")
elif shot == 184943:   # 50 kHz, DiMES experiment
    time_range = [2, 5]
    lowcut     = 1000
    highcut    = 2000
    sep_chord  = "FS4MIDDA"
    print("Warning: Tough sell...")
elif shot == 186827:   # 50 kHz, USN density ramp
    time_range = [2, 5]
    lowcut     = 1000
    highcut    = 2000
    sep_chord  = "FS5MIDDA"
    print("Warning: No real blob component in signal!")
elif shot == 186968:   # 50 kHz, LSN floor detachment experiment
    time_range = [2, 3.2]   # Density ramp starts right before 3.5s, can include if you want.
    lowcut     = 1550
    highcut    = 1900
    sep_chord  = "FS5MIDDA"


# Approximate locations of each midplane chord.
if mode == "midplane":
    midfs_rs = [2.208, 2.220, 2.233, 2.246, 2.258, 2.269, 2.281, 2.293]
    signals = ["FS{}MIDDA".format(i) for i in range(1,9)]
    midfs_rs = {s:r for s, r in zip(signals, midfs_rs)}
elif mode == "lower":
    signals = ["FS0{}".format(i) for i in range(0, 9)]
data = {}; good_signals = []; good_midfs_rs = {}
for signal in signals:
    print(signal)
    try:
        conn = MDSplus.Connection("atlas.gat.com")
        gaobj = gadata(signal, shot, connection=conn)
        t = gaobj.xdata / 1000
        y = gaobj.zdata
        mask = np.logical_and(t>time_range[0], t<time_range[1])
        t = t[mask]
        y = y[mask]

        # Fourier transform. Calculate average sample rate (should be ~constant).
        sample_rate = 1 / np.mean([t[i+1] - t[i] for i in range(1, len(t)-1)])
        print("  {:d} Hz".format(int(sample_rate)))
        yf = fft(y)
        xf = fftfreq(len(t), 1 / sample_rate)

        # Divide by the mean from 1000 Hz onwards as a way of normalizing.
        avg = np.abs(yf[xf > 1000]).mean()

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

        data[signal] = {"xf":xf, "yf":np.abs(yf) / avg, "rate":sample_rate,
            "x":t, "y":y, "yfilt":yfilt}
        good_signals.append(signal)
    except:
        print("  Skipped!")

if mode == "midplane":
    print("Don't trust anything below {:.2e} s".format(1 / (sample_rate * 0.5)))
    def lag_finder(y1, y2, sr):
        n = len(y1)

        corr = correlate(y2, y1, mode='same') / np.sqrt(correlate(y1, y1, mode='same')[int(n/2)] * correlate(y2, y2, mode='same')[int(n/2)])

        delay_arr = np.linspace(-0.5*n/sr, 0.5*n/sr, n)
        delay = delay_arr[np.argmax(corr)]
        #print('y2 is {:.3e} behind y1'.format(delay))

        #plt.figure()
        #plt.plot(delay_arr, corr)
        #plt.title('Lag: ' + str(np.round(delay, 3)) + ' s')
        #plt.xlabel('Lag')
        #plt.ylabel('Correlation coeff')
        #plt.show()

        return delay

    # Get the lags, for better or worse.
    print("Possible blob velocities (m/s)")
    for signal in good_signals:
        delay = lag_finder(data[sep_chord]["yfilt"], data[signal]["yfilt"], sample_rate)
        dist = midfs_rs[signal] -  midfs_rs[sep_chord]
        if delay == 0.0:
            print("{}: -----".format(signal))
        else:
            print("{}: {:.1f}".format(signal, dist / delay))
        data[signal]["lag"] = delay

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
