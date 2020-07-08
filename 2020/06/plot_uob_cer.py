import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import savgol_filter


# Load in all the data.
fname = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, " + \
        "Sheets and Documents/2020/06/startup_cer.xlsx"
df = pd.read_excel(fname, sheet_name="CER")

# Pull out data for the V chords.
v_r = df["R (cm).2"].values
#v_amp = df["Brightness (ph/s/m2/sr).2"].values
v_amp = df["CVI nz (m-3).2"].values
v_times = df["Time Range (ms).2"].values

# Indices for each time range.
idx_1000_1500 = np.where(v_times=='"1000-1500"')[0]
idx_1500_2500 = np.where(v_times=='"1500-2500"')[0]
idx_2500_3500 = np.where(v_times=='"2500-3500"')[0]
idx_3500_4500 = np.where(v_times=='"3500-4500"')[0]

# Sort each time subset accordiing to the R values.
idx_1000_1500 = [x for _,x in sorted(zip(v_r[idx_1000_1500],idx_1000_1500))]
idx_1500_2500 = [x for _,x in sorted(zip(v_r[idx_1500_2500],idx_1500_2500))]
idx_2500_3500 = [x for _,x in sorted(zip(v_r[idx_2500_3500],idx_2500_3500))]
idx_3500_4500 = [x for _,x in sorted(zip(v_r[idx_3500_4500],idx_3500_4500))]

# Plotting commands.
lw = 3
fig, ax = plt.subplots()
cmap = cm.get_cmap("Oranges")
ax.axvline(225.5, color="k", linestyle="--", lw=lw)
#ax.plot(v_r[idx_1000_1500], v_amp[idx_1000_1500], color=cmap(0.25), lw=lw, label="1000-1500")
ax.plot(v_r[idx_1500_2500], v_amp[idx_1500_2500], color=cmap(0.50), lw=lw, label="1500-2500")
ax.plot(v_r[idx_2500_3500], v_amp[idx_2500_3500], color=cmap(0.75), lw=lw, label="2500-3500")
ax.plot(v_r[idx_3500_4500], v_amp[idx_3500_4500], color=cmap(1.00), lw=lw, label="3500-4500")
ax.set_xlabel("R (cm)", fontsize=16)
#ax.set_ylabel("Amplitude (ph/s/m2/sr)", fontsize=16)
ax.set_ylabel("C6+ Density (m-3)", fontsize=16)
ax.legend(fontsize=12)
fig.tight_layout()
fig.show()
