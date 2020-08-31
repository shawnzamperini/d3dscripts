# Simple script to just print out the LAMS counts at each RBS location. The fits
# can then be done in Excel since it's easier to visualize and organize there.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Where we treat centerline as in the 2D LAMS scans.
pol = 2.5  # mm

# Just choose which probe we're doing here.
rbs_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Collector Probe Excel Sheets/B7.xlsx"
lams_pathD  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BD07_Map_Analysis.xlsx"
lams_pathU  = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Polodial_Scans/New Map Script Results/BU07_Map_Cal_Analysis.xlsx"

# Load the RBS and LAMS data.
rbs_df = pd.read_excel(rbs_path)
lamsD_df = pd.read_excel(lams_pathD, sheet_name="MapData")
lamsU_df = pd.read_excel(lams_pathU, sheet_name="MapData")

# Pull out the RBS locs and W data, removing nans.
def pull_rbs(df, side):
    """
    side: Either "U" or "D"
    """

    # Grab data.
    locs = df["Distance from Tip {} (cm)".format(side)].values
    w = df["W Areal Density {} (1e15 W/cm2)".format(side)].values

    # Filter nans.
    idx = ~np.isnan(locs)
    locs = locs[idx]
    w = w[idx]

    return locs, w

rbs_locsD, wD = pull_rbs(rbs_df, "D")
rbs_locsU, wU = pull_rbs(rbs_df, "U")

# Pull out the centerline LAMS data for each side.
idxD = lamsD_df["z Location [mm]"] == pol
lams_locsD = lamsD_df["Axial Location [mm]"][idxD].values / 10 # mm to cm
lams_totWD = lamsD_df["Total W"][idxD].values
idxU = lamsU_df["z Location [mm]"] == pol
lams_locsU = lamsU_df["Axial Location [mm]"][idxU].values / 10 # mm to cm
lams_totWU = lamsU_df["Total W"][idxU].values

# For each RBS location, find what the nearest point in LAMS is for counts.
calib_lamsD = np.zeros(len(rbs_locsD))
for i in range(0, len(rbs_locsD)):
    l = rbs_locsD[i]
    dist = np.abs(lams_locsD - l)
    close_idx = np.where(dist == dist.min())[0]
    calib_lamsD[i] = lams_totWD[close_idx]

calib_lamsU = np.zeros(len(rbs_locsU))
for i in range(0, len(rbs_locsU)):
    l = rbs_locsU[i]
    dist = np.abs(lams_locsU - l)
    close_idx = np.where(dist == dist.min())[0]
    calib_lamsU[i] = lams_totWU[close_idx]

# Save to a csv.
np.savetxt("lams_calib_tmp.csv", np.vstack((rbs_locsD, wD, calib_lamsD, rbs_locsU, wU, calib_lamsU)).T, delimiter=",")

# Plotting commands.
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
ax1b = ax1.twinx()
ax1.plot(rbs_locsD, wD, '.', ms=10, color="tab:cyan", mec="k")
ax1b.plot(lams_locsD, lams_totWD, color="tab:pink")
ax1.set_xlabel("Distance from tip (cm)", fontsize=16)
ax1b.set_ylabel("LAMS Counts", fontsize=16)
ax1.set_ylabel("W Areal Density D", fontsize=16)
ax2.plot(wD, calib_lamsD, '.', ms=10, color="tab:red")
ax2.set_xlabel("W Areal Density (1e15 w/cm2)", fontsize=16)
ax2.set_ylabel("LAMS Counts", fontsize=16)
ax3b = ax3.twinx()
ax3.plot(rbs_locsU, wU, '.', ms=10, color="tab:cyan", mec="k")
ax3b.plot(lams_locsU, lams_totWU, color="tab:pink")
ax3.set_xlabel("Distance from tip (cm)", fontsize=16)
ax3b.set_ylabel("LAMS Counts", fontsize=16)
ax3.set_ylabel("W Areal Density U", fontsize=16)
ax4.plot(wU, calib_lamsU, '.', ms=10, color="tab:red")
ax4.set_xlabel("W Areal Density (1e15 w/cm2)", fontsize=16)
ax4.set_ylabel("LAMS Counts", fontsize=16)
fig.tight_layout()
fig.show()
