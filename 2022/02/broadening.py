import numpy as np
from gadata import gadata
import MDSplus as mds
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


exclude = ["184170_1", "184170_2", "184173_2", "184174_2", "184177_2",
  "184178_2", "184179_2", "184176_1", "184182_2", "184264_1", "184264_2",
  "184266_2", "184268_1", "184268_2", "184270_1", "184270_2"]

# Some global variables.
nesep_path = "/Users/zamperini/My Drive/Research/Documents/2022/02/avg_neseps.xlsx"
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
conn = mds.Connection('atlas.gat.com')
nesep_df = pd.read_excel(nesep_path)

def exp_fit(x, a, b):
    return a * np.exp(-b * x)

rminrseps = {}
norm_nes = {}
fs_avgs = {}
machs = {}
dir = {}
avg_machs = {}
lambda_nes = {}
fitxs = {}
popts = {}

for i in range(0, len(nesep_df)):
    plunge = nesep_df["Shot_plunge"].iloc[i]
    print(plunge)
    if plunge in ["184270_1", "184270_2"]:
        print(" Skipped")
        continue
    shot = int(plunge.split("_")[0])
    start = nesep_df["Start"].iloc[i]
    end = nesep_df["End"].iloc[i]
    nesep = nesep_df["ne_sep"].iloc[i]

    rcp = pd.read_excel(rcp_path, sheet_name="MP{}".format(plunge))
    rcp_rminrsep = rcp["Dist. from sep."]
    rcp_ne = rcp["ne (1e18 m-3)"]

    rminrseps[plunge] = rcp_rminrsep
    norm_nes[plunge] = rcp_ne / nesep * 1e18
    dir[plunge] = nesep_df["Direction"].iloc[i]
    machs[plunge] = rcp["Mach"]

    # Average Mach number using only data within first 3 cm of separatrix.
    # Some pkunges are manually calculated.
    if plunge == "184175_1":
        avg_machs[plunge] = 0.101
    else:
        near = np.where(rcp_rminrsep<=0.03)[0]
        avg_machs[plunge] = rcp["Mach"][near].mean()

    # Fit to the ne data between 2-5 cm.
    if plunge == "184176_1":
        fit_start = 0.02
        fit_end = 0.05
    #if plunge in ["184154_1", "184154_2"]:
    #    fit_start = 0.02
    #    fit_end = 0.05
    else:
        fit_start = 0.02
        fit_end = 0.05
    fit_idx = np.where(np.logical_and(rcp_rminrsep >= fit_start, rcp_rminrsep <= fit_end))[0]
    fitx = rcp_rminrsep[fit_idx]
    fity = norm_nes[plunge][fit_idx]
    #try:
    popt, pcov = curve_fit(exp_fit, fitx, fity, maxfev=2000, bounds=[(0.0, 0.0), (np.inf, np.inf)], p0=(0.8, 10))
    lambda_ne = 1 / popt[1]
    #except:
    #print(" Error! Unable to fit exponential.")
    #lambda_ne = np.nan
    lambda_nes[plunge] = lambda_ne
    fitxs[plunge] = fitx
    popts[plunge] = popt

    # The ne / nesep value at

    fs02up = gadata("FS02UP", shot, connection=conn)
    fstime = fs02up.xdata
    fsdata = fs02up.zdata
    idx = np.where(np.logical_and(fstime>=start, fstime<=end))
    avg_fs = fsdata[idx].mean()

    fs_avgs[plunge] = avg_fs

max_for = 0
min_for = 1e40
max_rev = 0
min_rev = 1e40
for plunge, val in fs_avgs.items():
    if plunge in exclude:
        continue
    if dir[plunge] == "Forward":
        if val > max_for:
            max_for = val
        if val < min_for:
            min_for = val
    if dir[plunge] == "Reverse":
        if val > max_rev:
            max_rev = val
        if val < min_rev:
            min_rev = val

cmap = matplotlib.cm.get_cmap("inferno")
for_colors = {}
rev_colors = {}
for plunge, val in fs_avgs.items():
    if plunge in exclude:
        continue
    if dir[plunge] == "Forward":
        for_colors[plunge] = (val - min_for) / (max_for - min_for) * 0.8
    elif dir[plunge] == "Reverse":
        rev_colors[plunge] = (val - min_rev) / (max_rev - min_rev) * 0.8


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 5), sharex=True)

for plunge in rminrseps.keys():
    if plunge in exclude:
        continue
    if dir[plunge] == "Forward":
        ax1.plot(rminrseps[plunge], norm_nes[plunge], color=cmap(for_colors[plunge]),
          label=plunge)
        ax3.plot(rminrseps[plunge], machs[plunge], color=cmap(for_colors[plunge]))
    elif dir[plunge] == "Reverse":
        ax2.plot(rminrseps[plunge], norm_nes[plunge], color=cmap(rev_colors[plunge]),
          label=plunge)
        ax4.plot(rminrseps[plunge], machs[plunge], color=cmap(rev_colors[plunge]))

ax1.set_ylim([0.05, 1.5])
ax2.set_ylim([0.05, 1.5])
ax1.set_xlabel("R-Rsep (cm)")
ax2.set_xlabel("R-Rsep (cm)")
ax1.set_ylabel("ne / nesep")
ax1.set_yscale("log")
ax2.set_yscale("log")
#ax1.legend()
#ax2.legend()
fig.tight_layout()
fig.show()

# Reverse BT only for now.
mach_colors = {}
max_mach = 0
min_mach = 999
mach_cmap = matplotlib.cm.get_cmap("Reds")
for plunge, val in avg_machs.items():
    if plunge in exclude:
        continue
    if dir[plunge] == "Reverse":
        if val > max_mach:
            max_mach = val
        if val < min_mach:
            min_mach = val

# Since we are using a diverging colormap, set the max/min to the same (just minus).
#if abs(max_mach) > abs(min_mach):
#    min_mach = -max_mach
#else:
#    max_mach = -min_mach

mach_colors = {}
for plunge, val in avg_machs.items():
    if plunge in exclude:
        continue
    if dir[plunge] == "Reverse":
        mach_colors[plunge] = (val - min_mach) / (max_mach - min_mach)

fig, ax = plt.subplots()
x = []
y = []
c = []
z = []
p = []
for plunge in rminrseps.keys():
    if plunge in exclude:
        continue
    if dir[plunge] == "Forward":
        continue
    x.append(fs_avgs[plunge])
    y.append(lambda_nes[plunge])
    c.append(mach_cmap(mach_colors[plunge]))
    z.append(avg_machs[plunge])
    p.append(plunge)
    ax.annotate(plunge, (fs_avgs[plunge], lambda_nes[plunge]))
sc = ax.scatter(x, y, c=z, cmap="Reds", vmin=0, vmax=0.35, edgecolors="k")
ax.grid()
cbar = fig.colorbar(sc)
cbar.set_label("Mach Number between R-Rsep = 0-3 cm")
ax.set_xlabel("FS02UP")
ax.set_ylabel(r"$\mathdefault{\lambda_{ne}}$ between R-Rsep = 2-5 cm")
fig.tight_layout()
fig.show()

fig, axs = plt.subplots(5, 5, figsize=(12, 9))
axs = axs.flatten()

i = 0; j = 0
while i < len(axs):
#for i in range(0, len(fitxs)):
    print(i)
    if j >= len(fitxs.keys()):
        break
    plunge = list(fitxs.keys())[j]
    if plunge in exclude:
        j += 1
        continue
    if dir[plunge] == "Forward":
        j += 1
        continue
    axs[i].scatter(rminrseps[plunge], norm_nes[plunge], label=plunge, s=15)
    axs[i].plot(fitxs[plunge], exp_fit(fitxs[plunge], *popts[plunge]), color="k")
    axs[i].legend(loc="upper right", fontsize=10)
    axs[i].set_yscale("log")
    axs[i].set_xlim([0.0, 0.12])
    axs[i].set_ylim([0.05, 5])
    axs[i].grid()
    i += 1
    j += 1

fig.tight_layout()
fig.show()
