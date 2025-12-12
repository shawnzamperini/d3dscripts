import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from scipy.signal import medfilt


# Path to Excel file with data from flan_pygkyl of NT data. Need to mount 
# Google Drive with: sudo mount -t drvfs G: /mnt/g
path = "/mnt/g/My Drive/Research/Documents/2025/flan_nt_net_influx.xlsx"
df = pd.read_excel(path, skiprows=1, sheet_name="Sheet2")  # nt_v3_coll is here

# Pull out data
r = df["r (m)"].values
z = df["z (m)"].values
s = df["s (m)"].values
nw_coll = df["nw_coll"].values
Epol = df["Epol (V/m)"].values
Btor = df["Btor (T)"].values

# ExB drift across LCFS
ExB_lcfs = Epol / Btor

# Plot of W density
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
ax1.axvline(1.0, color="k", linestyle="--", lw=2)
ax1.axvline(3.06, color="k", linestyle="--", lw=2)
ax1.plot(s, nw_coll, color="k", lw=3)
ax1.plot(s, nw_coll, color="tab:red", lw=2)
#ax1.set_xlabel("Distance along separatrix (m)", fontsize=16)
ax1.set_ylabel("W Density (arb.)", fontsize=16)
ax2.axvline(1.0, color="k", linestyle="--", lw=2)
ax2.axvline(3.06, color="k", linestyle="--", lw=2)
ax2.plot(s, ExB_lcfs, color="k", lw=3)
ax2.plot(s, ExB_lcfs, color="tab:red", lw=2)
ax2.set_xlabel("Distance along separatrix (m)", fontsize=16)
ax2.set_ylabel("ExB_lcfs (m/s)", fontsize=16)
fig.tight_layout()
fig.show()

# We apply a median filter to the rate data to handle the noise
window = 21
rate_nocoll = medfilt(df["rate_nocoll"].values, window)
rate_coll = medfilt(df["rate_coll"].values, window)

# Determine min/max among both simulations for colorbar
min_rate = min(rate_nocoll.min(), rate_coll.min())
max_rate = max(rate_nocoll.max(), rate_coll.max())
min_ExB = ExB_lcfs.min()
max_ExB = ExB_lcfs.max()

# Then the vmin/vmax will span the range just so zero is in the middle of the
# colorbar in the plot
clip = 0.9
vmin = -max(abs(min_rate), abs(max_rate)) * clip
vmax = max(abs(min_rate), abs(max_rate)) * clip
vmin_ExB = -max(abs(min_ExB), abs(max_ExB)) * clip
vmax_ExB = max(abs(min_ExB), abs(max_ExB)) * clip

# Create segments
points = np.array([r, z]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Create LineCollection
cmap = "coolwarm"
lw = 16
norm = plt.Normalize(vmin, vmax)
lc_nocoll = LineCollection(segments, cmap=cmap, norm=norm)
lc_coll = LineCollection(segments, cmap=cmap, norm=norm)
lc_nocoll.set_array(rate_nocoll)
lc_coll.set_array(rate_coll)
lc_nocoll.set_linewidth(lw)
lc_coll.set_linewidth(lw)
lc_nocoll.set_antialiased(True)

# Again for ExB drift
norm_ExB = plt.Normalize(vmin_ExB, vmax_ExB)
lc_ExB = LineCollection(segments, cmap=cmap, norm=norm_ExB)
lc_ExB.set_array(ExB_lcfs)
lc_ExB.set_linewidth(lw)
lc_ExB.set_antialiased(True)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7), sharex=True, sharey=True)
ax1.set_aspect("equal")
ax2.set_aspect("equal")
ax1.set_title("Collisions OFF", fontsize=20)
ax2.set_title("Collisions ON", fontsize=20)
ax1.set_xlabel("R (m)", fontsize=18)
ax2.set_xlabel("R (m)", fontsize=18)
ax1.set_ylabel("Z (m)", fontsize=18)
ax1.tick_params(which="both", labelsize=16)
ax2.tick_params(which="both", labelsize=16)
ax1.plot(r, z, color="k", lw=lw+3)
ax1.plot(r, z, color="w", lw=lw)
ax2.plot(r, z, color="k", lw=lw+3)
ax2.plot(r, z, color="w", lw=lw)
ax1.add_collection(lc_nocoll)
ax2.add_collection(lc_coll)
ax1.autoscale()
ax2.autoscale()

# Adjust layout to make space for colorbar
fig.subplots_adjust(right=0.85)

# Add colorbar outside the figure
cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(lc_coll, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Rate of W Exiting Core (arb.)', fontsize=18)

fig.show()

# For the ExB drift
fig, ax1 = plt.subplots(figsize=(6, 7))
ax1.set_aspect("equal")
ax1.set_title("ExB Drift Inwards", fontsize=20)
ax1.set_xlabel("R (m)", fontsize=18)
ax1.set_ylabel("Z (m)", fontsize=18)
ax1.tick_params(which="both", labelsize=16)
ax1.plot(r, z, color="k", lw=lw+3)
ax1.plot(r, z, color="w", lw=lw)
ax1.add_collection(lc_ExB)
ax1.autoscale()

# Adjust layout to make space for colorbar
fig.subplots_adjust(right=0.85)

# Add colorbar outside the figure
cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(lc_ExB, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('ExB Drift Inward (m/s)', fontsize=18)

fig.show()
