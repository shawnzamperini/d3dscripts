import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import savgol_filter


# Load in all the data.
fname = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, " + \
        "Sheets and Documents/2020/06/startup_cer.xlsx"
df = pd.read_excel(fname, sheet_name="TS", skiprows=1)

# Plotting commands.
lw = 3
fig, ax = plt.subplots()
cmap = cm.get_cmap("Oranges")
ax.axvline(0.725, color="k", lw=lw, linestyle="--")
#ax.plot(df["Z (m)"], df["Te (eV)"], label="1000-1500", color=cmap(0.25))
ax.plot(df["Z (m)"], df["Te (eV).1"], label="1500-2500", color=cmap(0.50))
ax.plot(df["Z (m)"], df["Te (eV).2"], label="2500-3500", color=cmap(0.75))
ax.plot(df["Z (m)"], df["Te (eV).3"], label="3500-4500", color=cmap(1.00))
ax.legend(fontsize=12)
ax.set_xlabel("Z (m)", fontsize=16)
ax.set_ylabel("Te (eV)", fontsize=16)
fig.tight_layout()
fig.show()
