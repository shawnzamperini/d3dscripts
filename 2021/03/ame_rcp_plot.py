import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import savgol_filter


# Paths to the plunge data for each RCP.
mp_fav_dfs = []; mp_unf_dfs = []
path1 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_1.tab".format("MP", 184267)
path2 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_2.tab".format("MP", 184267)
df1 = pd.read_csv(path1, sep="\t")
df2 = pd.read_csv(path2, sep="\t")
mp_fav_dfs.append(df1)
mp_fav_dfs.append(df2)
path1 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_1.tab".format("MP", 184527)
path2 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_2.tab".format("MP", 184527)
df1 = pd.read_csv(path1, sep="\t")
df2 = pd.read_csv(path2, sep="\t")
mp_unf_dfs.append(df1)
mp_unf_dfs.append(df2)

xp_fav_dfs = []; xp_unf_dfs = []
path1 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_1.tab".format("XP", 184267)
path2 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_2.tab".format("XP", 184267)
df1 = pd.read_csv(path1, sep="\t")
df2 = pd.read_csv(path2, sep="\t")
xp_fav_dfs.append(df1)
xp_fav_dfs.append(df2)
path1 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_1.tab".format("XP", 184527)
path2 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_2.tab".format("XP", 184527)
df1 = pd.read_csv(path1, sep="\t")
df2 = pd.read_csv(path2, sep="\t")
xp_unf_dfs.append(df1)
xp_unf_dfs.append(df2)

lw = 3

fig, ax1 = plt.subplots()
ax1.plot(mp_fav_dfs[0]["Rho"], mp_fav_dfs[0]["Machn"], color="tab:red", lw=lw)
ax1.plot(mp_fav_dfs[1]["Rho"], mp_fav_dfs[1]["Machn"], color="tab:red", lw=lw)
ax1.plot(mp_unf_dfs[0]["Rho"], mp_unf_dfs[0]["Machn"], color="tab:purple", lw=lw)
ax1.plot(mp_unf_dfs[1]["Rho"], mp_unf_dfs[1]["Machn"], color="tab:purple", lw=lw)
ax1.set_title("MiMES RCP", fontsize=16)
ax1.set_xlabel("Rho", fontsize=16)
ax1.set_ylabel("Mach Number", fontsize=16)
ax1.grid()
ax1.set_ylim([-1, 1])
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
custom_lines = [Line2D([0], [0], color="tab:red"),
                Line2D([0], [0], color="tab:purple")]
ax1.legend(custom_lines, ["Favorable", "Unfavorable"], fontsize=12)
ax1.axhline(0.0, color="k", linestyle="-", lw=4)
fig.tight_layout()
fig.show()

fig, ax1 = plt.subplots()
ax1.plot(xp_fav_dfs[0]["Rho"][9:-9], xp_fav_dfs[0]["Machn"][9:-9], color="tab:red", lw=lw, alpha=1.0)
ax1.plot(xp_fav_dfs[1]["Rho"][9:-9], xp_fav_dfs[1]["Machn"][9:-9], color="tab:red", lw=lw, alpha=1.0)
#ax1.plot(xp_fav_dfs[0]["Rho"][9:-9], savgol_filter(xp_fav_dfs[0]["Machn"][9:-9], 41, 2), color="tab:red", lw=lw)
#ax1.plot(xp_fav_dfs[1]["Rho"][9:-9], savgol_filter(xp_fav_dfs[1]["Machn"][9:-9], 41, 2), color="tab:red", lw=lw)
ax1.plot(xp_unf_dfs[0]["Rho"][5:-3], xp_unf_dfs[0]["Machn"][5:-3], color="tab:purple", lw=lw)
ax1.plot(xp_unf_dfs[1]["Rho"][5:-3], xp_unf_dfs[1]["Machn"][5:-3], color="tab:purple", lw=lw)
ax1.set_title("X-point RCP", fontsize=16)
ax1.set_xlabel("Rho", fontsize=16)
ax1.set_ylabel("Mach Number", fontsize=16)
ax1.grid()
ax1.set_ylim([-1, 1])
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
custom_lines = [Line2D([0], [0], color="tab:red"),
                Line2D([0], [0], color="tab:purple")]
ax1.legend(custom_lines, ["Favorable", "Unfavorable"], fontsize=12)
ax1.axhline(0.0, color="k", linestyle="-", lw=4)
fig.tight_layout()
fig.show()
