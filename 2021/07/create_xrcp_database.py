import MDSplus
import pandas as pd
import numpy  as np
from   gadata import gadata
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib as mpl


# Some constants.
rho_limit = 1.055  # There's a clear dropoff in Te signal after this point.
unf_shots = np.append(np.arange(184151, 184164), np.arange(184524, 184541))
fav_shots = np.append(np.arange(184164, 184185), np.arange(184262, 184275))

# Load in the Excel file with all the RCP data, and setup an MDSplus Connection
# object.
rcp_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"
dfs = pd.read_excel(rcp_path, sheet_name=None)
conn = MDSplus.Connection("localhost")

# Helper function to get time-average values for a shot.
def time_avg_value(shot, tag, time_start, time_end):

    # Special use case where we need to interpolate onto a common time basis.
    if tag == "PSOL":
        pinj = gadata("PINJ", shot, connection=conn)
        prad_core = gadata("PRAD_CORE", shot, connection=conn)
        f_pinj = interp1d(pinj.xdata, pinj.zdata, bounds_error=False, fill_value=0)
        psol = f_pinj(prad_core.xdata) * 1000 - prad_core.zdata
        keep = np.logical_and(prad_core.xdata>=time_start, prad_core.xdata<=time_end)
        avg = psol[keep].mean()

    else:
        gaobj = gadata(tag, shot, connection=conn)
        keep = np.logical_and(gaobj.xdata>=time_start, gaobj.xdata<=time_end)
        avg = gaobj.zdata[keep].mean()

    return avg

fig, ax = plt.subplots()

# Go through one plunge at a time.
fav_count = 0
unf_count = 0
all_df = pd.DataFrame()
for rcp in dfs.keys():

    # Only do XP probes.
    if rcp[:2] != "XP":
        continue

    # Can pull out the shot, as well at the time ranges.
    rcp_df = dfs[rcp]
    shot = int(rcp.split("P")[1].split("_")[0])
    time_start = rcp_df["Time(ms)"].min()
    time_end   = rcp_df["Time(ms)"].max()
    print(shot)

    if shot in fav_shots:
        fav_count += 1
        rcp_df["direction"] = np.full(len(rcp_df), "favorable")
    elif shot in unf_shots:
        unf_count += 1
        rcp_df["direction"] = np.full(len(rcp_df), "unfavorable")
    else:
        print("Error: {} not identified as either direction.".format(shot))

    # Mask the data beyond rho_limit.
    mask = rcp_df["Rho"] < rho_limit

    # Add average plasma values during the plunge.
    rcp_df["ne_ped"] = np.full(len(rcp_df), time_avg_value(shot, "PRMTAN_NEPED", time_start, time_end))
    rcp_df["te_ped"] = np.full(len(rcp_df), time_avg_value(shot, "PRMTAN_TEPED", time_start, time_end))
    rcp_df["pinj"]   = np.full(len(rcp_df), time_avg_value(shot, "PINJ", time_start, time_end))
    rcp_df["tinj"]   = np.full(len(rcp_df), time_avg_value(shot, "TINJ", time_start, time_end))
    rcp_df["densv2"] = np.full(len(rcp_df), time_avg_value(shot, "DENSV2", time_start, time_end))
    rcp_df["psol"]   = np.full(len(rcp_df), time_avg_value(shot, "PSOL", time_start, time_end))

    # Assign identifier and put into master DataFrame.
    rcp_df["id"] = np.full(len(rcp_df), rcp)
    all_df = all_df.append(rcp_df[mask])

    ax.plot(rcp_df["Rho"], rcp_df["Machn"])

print("Favorable:   {}".format(fav_count))
print("Unfavorable: {}".format(unf_count))

ax.set_xlabel("Rho")
ax.set_ylabel("Mach")
fig.tight_layout()
fig.show()

# Save as an Excel.
all_df.to_excel("xrcp_data.xlsx")

# Little plot down here to keep things in one place.
fig, ax = plt.subplots()

color_boss = "pinj"
all_df["color_scale"] = (all_df[color_boss] - all_df[color_boss].min()) / (all_df[color_boss].max() - all_df[color_boss].min())
cmap = plt.get_cmap('Reds')

for id in all_df["id"].unique():
    df = all_df[all_df["id"] == id]
    color = cmap(df["color_scale"].mean())
    ax.plot(df["Rho"], df["Machn"], color=color)

    # Additional lines to identify a couple lines.
    if id == "XP184537_1":
        ax.plot(df["Rho"], df["Machn"], color="g")
    if id == "XP184270_2":
        ax.plot(df["Rho"], df["Machn"], color="b")

ax.set_xlabel("Rho")
ax.set_ylabel("Mach")
fig.tight_layout()
fig.show()
