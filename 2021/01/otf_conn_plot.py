import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


shots = [167530, 167534, 167536, 167279, 167320, 167321, 167322, 167353, 167405, 167355, 167358, 167377, 167380]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12,5))

for shot in shots:

    # Load the connection length data from the Excel sheet.
    path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/Connection Lengths/{}/{}.xlsx".format(shot, shot)
    df = pd.read_excel(path, sheet_name="MAFOT OTF", skiprows=2)
    df2 = pd.read_excel(path, sheet_name="MAFOT ITF", skiprows=2)

    x = df["R-Rsep OMP (cm)"].values
    y = df["Connection Length (km)"].values * 1000  # km to m
    x2 = df2["R-Rsep OMP (cm)"].values
    y2 = df2["Connection Length (km)"].values * 1000  # km to m

    ax.plot(x, y, label=shot)
    ax2.plot(x2, y2, label=shot)

    # Also just print out the average connection length from 6-13.
    mask = np.logical_and(x>=6, x<=13)
    avg_conn = np.mean(y[mask])
    print("{}: {:.2f}".format(shot, avg_conn))

ax.set_xlabel("R-Rsep OMP (cm)", fontsize=16)
ax2.set_xlabel("R-Rsep OMP (cm)", fontsize=16)
ax.set_ylabel("Distance to nearest limiter (m)", fontsize=16)
ax.legend(ncol=2)
ax.set_xlim([0, 15])
ax2.set_xlim([0, 15])
ax.set_ylim([1, 100])
ax2.set_ylim([1, 100])
ax.set_yscale("log")
ax2.set_yscale("log")
ax.tick_params(which='both', labelsize=12)
ax2.tick_params(which='both', labelsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
