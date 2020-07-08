import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# Load in all the data.
fname = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, " + \
        "Sheets and Documents/2020/06/startup_cer.xlsx"

df = pd.read_excel(fname, sheet_name="MDS")

uobx = df["UOB Time (s)"].values
uoby = df["UOB CAL"].values
mdsx = df["Time (s)"].values
#mdsy = df["Peak Intensity"].values
mdsy = df["L4 Intensity"].values

# Smoothed MDS data.
mdsy_smooth = savgol_filter(mdsy, 15, 3)

# Plotting commands.
fig, ax1 = plt.subplots()
ax1.plot(mdsx, mdsy, color='b', lw=3, alpha=0.25)
#ax1.plot(mdsx, mdsy_smooth, color='b', lw=3, label="MDS U5")
ax1.plot(mdsx, mdsy_smooth, color='b', lw=3, label="MDS L4")
ax2 = ax1.twinx()
ax2.plot(uobx, uoby, color='r', alpha=0.75, label="UOB")
ax2.set_xlim([0, 6])
ax1.set_zorder(10)
ax1.patch.set_visible(False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax1.set_xlabel("Time (s)", fontsize=16)
#ax1.set_ylabel("MDS U5 Intensity", fontsize=16)
ax1.set_ylabel("MDS L4 Intensity", fontsize=16)
#ax1.legend()
#ax2.legend()
fig.tight_layout()
fig.show()
