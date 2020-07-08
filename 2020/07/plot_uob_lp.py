import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import savgol_filter


fname = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, " + \
        "Sheets and Documents/2020/06/startup_cer.xlsx"
df = pd.read_excel(fname, sheet_name="LP")

wind_len = 35
fontsize = 16
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
for t in [1500, 2500, 3500]:
    ax1.axvline(t, color='k', linestyle="--", lw=3)
    ax2.axvline(t, color='k', linestyle="--", lw=3)
ax1.plot(df['17 Time'], df['17 Te'], color='tab:pink', alpha=0.25)
ax1.plot(df['19 Time'], df['19 Te'], color='tab:cyan', alpha=0.25)
ax1.plot(df['17 Time'], savgol_filter(df['17 Te'], wind_len, 3), color='tab:pink')
ax1.plot(df['19 Time'], savgol_filter(df['19 Te'], wind_len, 3), color='tab:cyan')
ax2.plot(df['17 Time'], df['17 jsat'], color='tab:pink', alpha=0.25)
ax2.plot(df['19 Time'], df['19 jsat'], color='tab:cyan', alpha=0.25)
ax2.plot(df['17 Time'], savgol_filter(df['17 jsat'], wind_len, 3), color='tab:pink', label="LP17")
ax2.plot(df['19 Time'], savgol_filter(df['19 jsat'], wind_len, 3), color='tab:cyan', label="LP19")
ax2.legend(loc="upper right", fontsize=12)
ax1.set_ylabel("Te (eV)", fontsize=fontsize)
ax2.set_ylabel("jsat (A/cm2)", fontsize=fontsize)
ax1.set_xlabel("Time (ms)", fontsize=fontsize)
ax2.set_xlabel("Time (ms)", fontsize=fontsize)
ax1.set_xlim([1000,5000])
ax1.set_ylim([0, 50])
ax2.set_xlim([1000,5000])
fig.tight_layout()
fig.show()
