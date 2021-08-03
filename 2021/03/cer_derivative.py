import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from gadata import gadata
from scipy.interpolate import interp1d


# We know what we're doing here so ignore these.
pd.options.mode.chained_assignment = None

# Load the csv file from OMFIT FIT data.
shot = 184271
path = "/mnt/c/Users/Shawn/Documents/d3d_work/Methane Experiment/n_12C6_{}.csv".format(shot)
df = pd.read_csv(path)

if shot == 184527:
    puff_start = 2800
    xlim       = [None, None]
    ylim       = [None, None]
    window     = 1501
    beam_time = [1400, 5400]
elif shot == 184535:
    puff_start = 2200
    xlim       = [None, None]
    ylim       = [None, None]
    window     = 2501
    beam_time  = [1400, 5400]
elif shot == 184271:
    puff_start = 2200
    xlim       = [1500, None]
    #xlim      = [None, None]
    #ylim      = [None, 3e14]
    ylim       = [-7e-6, 7e-6]
    window     = 2501
    beam_time  = [1400, 5400]
elif shot == 184272:
    puff_start = 2200
    xlim       = [None, None]
    ylim       = [None, None]
    window     = 1501
    beam_time  = [1400, 5400]

# Convert the string uncertainty data to floats.
df["n_12C6"] = df["n_12C6"].apply(lambda x: float(re.sub("[+].*[)]", "", x)[1:]))

# Cut of the earliest (meaningless) times so interpolation with ne works.
df = df[df["time"] > 150]

# Restrict to just data when the beams were on.
df = df[df["time"] > beam_time[0]]
df = df[df["time"] < beam_time[1]]

# Load density data and interpolate on the CER time scale.
ga = gadata("DENSV2", shot)
ne_time = ga.xdata
ne = ga.zdata
f = interp1d(ne_time, ne)
ne_int = f(df["time"])
df["n_12C6_ne"] = df["n_12C6"] / ne_int

# New column of the time derivative.
df["n_12C6_smooth"] = np.zeros(len(df))
df["n_12C6_smooth_ne"] = np.zeros(len(df))
df["d_n_12C6_dt"] = np.zeros(len(df))
df["d_n_12C6_dt_smooth"] = np.zeros(len(df))
for psin in df["psi_n"].unique():
    tmp_df = df[df["psi_n"] == psin]
    tmp_df["n_12C6_smooth_ne"] = savgol_filter(tmp_df["n_12C6_ne"], window, 2, mode="nearest")
    tmp_df["d_n_12C6_dt_smooth"] = np.gradient(tmp_df["n_12C6_smooth_ne"].values, tmp_df["time"].values)
    df.update(tmp_df)
#df["n_12C6_smooth_ne"] = df["n_12C6_smooth"] / ne_int

# Plotting.
fontsize = 16
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
ax2.axhline(0.0, color="k", linestyle="-")
ax1.axvline(puff_start, color="k", linestyle="--")
ax2.axvline(puff_start, color="k", linestyle="--")

count = 0
for psin in [0.1, 0.5, 0.8, 1.0, 1.1]:
    color = "C{}".format(count)
    tmp_df = df[df["psi_n"] == psin]
    ax1.plot(tmp_df["time"], tmp_df["n_12C6_ne"]*100, label=psin, alpha=0.5, color=color)
    ax1.plot(tmp_df["time"], tmp_df["n_12C6_smooth_ne"]*100, label=psin, color=color)
    ax2.plot(tmp_df["time"], tmp_df["d_n_12C6_dt_smooth"], label=psin, color=color)
    count += 1

ax1.set_xlim(xlim)
ax2.set_xlim(xlim)
ax1.set_ylim(0, 3)
ax2.set_ylim(ylim)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.legend(fontsize=12, loc="lower right")
ax2.grid()
ax1.set_xlabel("Time (ms)", fontsize=fontsize)
ax1.set_ylabel(r"$\mathdefault{f_C=n_C/\bar{n}_e} (\%)$", fontsize=fontsize)
ax2.set_xlabel("Time (ms)", fontsize=fontsize)
ax2.set_ylabel(r"$\mathdefault{df_C/dt}$", fontsize=fontsize)
fig.suptitle(shot, fontsize=20)
fig.tight_layout()
fig.show()
