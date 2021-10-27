import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

mD = (2 * 931.494 * 10**6) / (3*10**8)**2  # eV * s2 / m2

# All the RCP data.
xl_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/" + \
  "rcp_master_detailed.xlsx"
master_xl = pd.ExcelFile(xl_path)

# Before starting, let's go through and get a list of every shot and just
# assign color numbers to them so they can be used to match XRCP and MRCP lines.
shot_colors = {}
count = 0
for sheet_name in master_xl.sheet_names:
    shot = int(sheet_name[2:8])
    if shot not in shot_colors.keys():
        shot_colors[shot] = "C{}".format(count)
        count += 1

        # Only use 10 colors. So yeah some will be repeated.
        if count > 9:
            count = 0

# Indicate which shots were during which day.
day1 = np.arange(184151, 184164, dtype=int)  # Forward, Unfavorable
day2 = np.arange(184164, 184185, dtype=int)  # Reverse, Favorable
day3 = np.arange(184262, 184275, dtype=int)  # Reverse, Favorable
day4 = np.arange(184524, 184541, dtype=int)  # Forward, Unfavorable
day5 = np.arange(187098, 187124, dtype=int)  # Forward, Unfavorable

exclude_shots = [184176, 184182, 184183, 184184, 184264, 184266, 184268, 184270]

all_day2 = {"psin":[], "z":[], "vf":[], "Bt":[], "Te":[]}
all_day3 = {"psin":[], "z":[], "vf":[], "Bt":[], "Te":[]}
all_day4 = {"psin":[], "z":[], "vf":[], "Bt":[], "Te":[]}

for sheet_name in master_xl.sheet_names:
        shot = int(sheet_name[2:8])
        ptype = sheet_name[:2]

        if shot in exclude_shots:
            continue

        if ptype == "MP":
            continue

        rcp_df = master_xl.parse(sheet_name)
        psin = rcp_df["Psin"].values
        z = rcp_df["Z (cm)"].values
        vf = rcp_df["Vf"].values
        bt = rcp_df["Bt"].values
        te = rcp_df["Te (eV)"].values

        if shot in day2:
            all_day2["psin"] = np.append(all_day2["psin"], psin)
            all_day2["z"] = np.append(all_day2["z"], z/100)
            all_day2["vf"] = np.append(all_day2["vf"], vf)
            all_day2["Bt"] = np.append(all_day2["Bt"], bt)
            all_day2["Te"] = np.append(all_day2["Te"], te)
        elif shot in day3:
            all_day3["psin"] = np.append(all_day3["psin"], psin)
            all_day3["z"] = np.append(all_day3["z"], z/100)
            all_day3["vf"] = np.append(all_day3["vf"], vf)
            all_day3["Bt"] = np.append(all_day3["Bt"], bt)
            all_day3["Te"] = np.append(all_day3["Te"], te)
        elif shot in day4:
            all_day4["psin"] = np.append(all_day4["psin"], psin)
            all_day4["z"] = np.append(all_day4["z"], z/100)
            all_day4["vf"] = np.append(all_day4["vf"], vf)
            all_day4["Bt"] = np.append(all_day4["Bt"], bt)
            all_day4["Te"] = np.append(all_day4["Te"], te)

def exp_fit(x, a, b, c):
    return a * np.exp(b * x) + c

# Sort the data according to psin (z would be fine too I guess).
for dic in [all_day2, all_day3, all_day4]:
    sort_idx = np.argsort(dic["z"])
    dic["psin"] = np.array(dic["psin"])[sort_idx]
    dic["z"] = np.array(dic["z"])[sort_idx]
    dic["vf"] = np.array(dic["vf"])[sort_idx]
    dic["Bt"] = np.array(dic["Bt"])[sort_idx]
    dic["Te"] = np.array(dic["Te"])[sort_idx]

    # Smooth the Vf data some.
    dic["vf_sg"] = savgol_filter(dic["vf"], 101, 2)

    # Exponential fit.
    popt, pcov = curve_fit(exp_fit, dic["z"], dic["vf"], maxfev=5000)
    dic["vf_exp"] = exp_fit(dic["z"], *popt)

    # E = -dVf/dz
    #dic["E"] = -np.gradient(dic["vf_sg"], dic["z"])
    dic["E"] = -popt[1] * exp_fit(dic["z"], *popt)

    # Pol ExB velocity.
    dic["ExB_pol"] = dic["E"] * dic["Bt"]

    # Represent as a Mach number.
    dic["cs"] = np.sqrt((2 * dic["Te"]) / mD)
    dic["cs_sg"] = savgol_filter(dic["cs"], 101, 2)
    dic["ExB_pol_M"] = dic["ExB_pol"] / dic["cs_sg"]

fig, axs = plt.subplots(3, 3, figsize=(10,7), sharex=True)
axs = axs.flatten()

axs[0].set_xticks([-1.25, -1.20, -1.15, -1.10, -1.05])

for ax in axs:
    ax.grid()

axs[0].plot(all_day2["z"], all_day2["vf"])
axs[3].plot(all_day3["z"], all_day3["vf"])
axs[6].plot(all_day4["z"], all_day4["vf"])
axs[0].plot(all_day2["z"], all_day2["vf_exp"])
axs[3].plot(all_day3["z"], all_day3["vf_exp"])
axs[6].plot(all_day4["z"], all_day4["vf_exp"])

axs[1].plot(all_day2["z"], all_day2["E"])
axs[4].plot(all_day3["z"], all_day3["E"])
axs[7].plot(all_day4["z"], all_day4["E"])

axs[2].plot(all_day2["z"], all_day2["ExB_pol_M"])
axs[5].plot(all_day3["z"], all_day3["ExB_pol_M"])
axs[8].plot(all_day4["z"], all_day4["ExB_pol_M"])
#axs[2].plot(all_day2["z"], all_day2["cs_sg"])
#axs[5].plot(all_day3["z"], all_day3["cs_sg"])
#axs[8].plot(all_day4["z"], all_day4["cs_sg"])

axs[0].set_ylim([-90, 10])
axs[3].set_ylim([-90, 10])
axs[6].set_ylim([-90, 10])
axs[1].set_ylim([0, 2400])
axs[4].set_ylim([0, 2400])
axs[7].set_ylim([0, 2400])
axs[2].set_ylim([-0.1, 0.1])
axs[5].set_ylim([-0.1, 0.1])
axs[8].set_ylim([-0.1, 0.1])

axs[6].set_xlabel("Z (m)")
axs[7].set_xlabel("Z (m)")
axs[8].set_xlabel("Z (m)")
axs[0].set_title("Vf (V)")
axs[1].set_title("Er (V/m)")
axs[2].set_title("Poloidal ExB (m/s)")
axs[0].set_ylabel("Favorable")
axs[3].set_ylabel("Favorable")
axs[6].set_ylabel("Unfavorable")

fig.tight_layout()
fig.show()
