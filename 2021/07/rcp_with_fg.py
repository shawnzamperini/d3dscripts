import MDSplus
import pandas as pd
import numpy  as np
from   gadata import gadata
import matplotlib.pyplot as plt

rcp_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"
dfs = pd.read_excel(rcp_path, sheet_name=None)

unf_shots = np.append(np.arange(184151, 184164), np.arange(184524, 184541))
fav_shots = np.append(np.arange(184164, 184185), np.append(np.arange(184262, 184275), np.arange(187103, 187112)))

conn = MDSplus.Connection("localhost")

dens = {}; fg = {}; avg_machs = []; std_machs = []; unf = []; xp = []
for rcp in dfs.keys():
    shot = int(rcp.split("P")[1].split("_")[0])
    time_start = dfs[rcp]["Time(ms)"].min()
    time_end   = dfs[rcp]["Time(ms)"].max()

    # Get average density during plunge.
    gaobj = gadata("DENSV2", shot, connection=conn)
    #gaobj = gadata("PRMTAN_NEPED", shot, connection=conn)
    keep = np.logical_and(gaobj.xdata>=time_start, gaobj.xdata<=time_end)
    avg_ne = gaobj.zdata[keep].mean()

    # Get average IP during plunge.
    gaobj = gadata("IP", shot, connection=conn)
    keep = np.logical_and(gaobj.xdata.value>=time_start, gaobj.xdata.value<=time_end)
    avg_ip = gaobj.zdata.value[keep].mean()

    # Get everage minor radius.
    gaobj = gadata("AMINOR", shot, connection=conn)
    keep = np.logical_and(gaobj.xdata>=time_start, gaobj.xdata<=time_end)
    avg_a = gaobj.zdata[keep].mean()

    # Calculate average greenwald fraction.
    avg_ng = avg_ip * 10**(-6) / (np.pi * avg_a**2) * 1e20
    avg_fg = avg_ne / avg_ng
    fg[rcp] = avg_fg

min_fg = 99999; max_fg = 0
for rcp in dfs.keys():
    if fg[rcp] < min_fg:
        min_fg = fg[rcp]
    if fg[rcp] > max_fg:
        max_fg = fg[rcp]

plot_rho = True
cmap = plt.get_cmap('magma')
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharey=True, figsize=(8, 6))
for rcp in dfs.keys():

    shot = int(rcp[2:8])
    rho = dfs[rcp]["Rho"]
    mach = dfs[rcp]["Machn"]
    #mach = dfs[rcp]["Ne (E18 m-3)"]

    color_val = (fg[rcp] - min_fg) / (max_fg - min_fg)
    if rcp[:2] == "MP":
        r = dfs[rcp]["R(cm)"]
        if shot in unf_shots:
            if plot_rho:
                ax1.plot(rho, mach, color=cmap(color_val))
            else:
                ax1.plot(r, mach, color=cmap(color_val))
        elif shot in fav_shots:
            if plot_rho:
                ax3.plot(rho, mach, color=cmap(color_val))
            else:
                ax3.plot(r, mach, color=cmap(color_val))
        else:
            print("Error: {}".format(rcp))
    elif rcp[:2] == "XP":
        z = dfs[rcp]["Z(cm)"]
        if shot in unf_shots:
            if plot_rho:
                ax2.plot(rho, mach, color=cmap(color_val))
            else:
                ax2.plot(z, mach, color=cmap(color_val))
        elif shot in fav_shots:
            if plot_rho:
                ax4.plot(rho, mach, color=cmap(color_val))
            else:
                ax4.plot(z, mach, color=cmap(color_val))

ax1.set_title("Midplane")
ax2.set_title("X-Point")
ax3.set_xlabel("R (cm)")
ax4.set_xlabel("Z (cm)")
ax1.set_ylabel("Unfavorable Mach")
ax3.set_ylabel("Favorable Mach")
ax1.set_ylim([-1, 1])
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
fig1.tight_layout()
fig1.show()
