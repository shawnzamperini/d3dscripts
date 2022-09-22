from gadata import gadata
import numpy as np
import MDSplus
from scipy.signal import find_peaks, savgol_filter
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib as mpl


shots = [190160, 190161, 190162, 190163, 190167, 190169, 190170, 190171, 190173,
    190175, 190176, 190177, 190178, 190180, 190181, 190182, 190183, 190184]

signal = "FS02UP"
times = [2000, 5000]
width = 200

conn = MDSplus.Connection("atlas.gat.com")

results = {}
for shot in tqdm(shots):
    gaobj = gadata(signal, shot, connection=conn)
    time = gaobj.xdata
    fs = gaobj.zdata
    mask = np.logical_and(time>times[0], time<times[1])
    time = time[mask]
    fs = fs[mask]

    peaks = find_peaks(fs, distance=200, height=2e15)
    peak_idx = peaks[0]

    # Count number of peaks in each bin.
    bins = np.arange(times[0], times[1], width)
    counts = np.zeros(bins.shape)
    for i in range(0, len(peak_idx)):
        t = time[peak_idx[i]]
        for j in range(0, len(bins)-1):
            if t >= bins[j] and t < bins[j+1]:
                counts[j] += 1
    freqs = counts / (width/1000)
    bin_mids = [bins[i] + (bins[i+1] - bins[i])/2 for i in range(0, len(bins)-1)]

    # Load SPRED FE23 signal.
    gaobj2 = gadata("SPRED_FE23", shot, connection=conn)
    time2 = gaobj2.xdata
    fe23 = savgol_filter(gaobj2.zdata, 31, 2)

    # Grab nearest SPRED data at the center of each ELM bin.
    fe23_mids = []
    for mid_loc in bin_mids:
        idx = np.argmin(np.abs(time2-mid_loc))
        fe23_mids.append(fe23[idx])

    # Pedestal temprature as proxy for core temperature.
    try:
        gaobj3 = gadata("ECHPWRC", shot, connection=conn)
        time3 = gaobj3.xdata
        teped = savgol_filter(gaobj3.zdata, 31, 2)

        # Grab nearest SPRED data at the center of each ELM bin.
        teped_mids = []
        for mid_loc in bin_mids:
            idx = np.argmin(np.abs(time3-mid_loc))
            teped_mids.append(teped[idx])
    except:
        teped_mids = []

    results[shot] = {"freqs":freqs, "fe23_mids":fe23_mids, "teped_mids":teped_mids}

# Create colorscale from in to max of pedestal Te.
all_tepeds = []
for shot in results.keys():
    tepeds = results[shot]["teped_mids"]
    all_tepeds.append(tepeds)
all_tepeds = np.array(all_tepeds).flatten()
norm = mpl.colors.Normalize(vmin=all_tepeds.min(), vmax=all_tepeds.max())

# Plot it all.
fig, ax1 = plt.subplots(figsize=(5, 4))
for shot in results.keys():
    x = results[shot]["freqs"][:-1]
    y = results[shot]["fe23_mids"]
    z = results[shot]["teped_mids"]
    if z == []:
        continue
    ax1.scatter(x, y, s=20, zorder=10, c=norm(z), alpha=0.7)
ax1.grid(zorder=1)
ax1.set_xlim([0, 150])
ax1.set_ylim([2e11, 4e13])
ax1.set_yscale("log")
ax1.set_xlabel("ELM Frequency")
ax1.set_ylabel("SPRED_Fe23")
fig.tight_layout()
fig.show()
