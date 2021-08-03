import MDSplus
import pandas as pd
import numpy  as np
from   gadata import gadata
import matplotlib.pyplot as plt

rcp_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"
dfs = pd.read_excel(rcp_path, sheet_name=None)

conn = MDSplus.Connection("localhost")

dens = {}; fg = {}; avg_machs = []; std_machs = []; unf = []; xp = []
for rcp in dfs.keys():
    shot = int(rcp.split("P")[1].split("_")[0])
    time_start = dfs[rcp]["Time(ms)"].min()
    time_end   = dfs[rcp]["Time(ms)"].max()

    # Get average density during plunge.
    #gaobj = gadata("DENSV2", shot, connection=conn)
    gaobj = gadata("PRMTAN_NEPED", shot, connection=conn)
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

    # Average Mach number.
    avg_machs.append(dfs[rcp]["Machn"].mean())
    #avg_machs.append(dfs[rcp]["Machn"].median())
    std_machs.append(dfs[rcp]["Machn"].std())

    # Get the BT direction.
    gaobj = gadata("BT", shot, connection=conn)
    keep = np.logical_and(gaobj.xdata>=time_start, gaobj.xdata<=time_end)
    avg_bt = gaobj.zdata[keep].mean()
    if avg_bt > 0:
        unf.append(False)
    else:
        unf.append(True)

    # Assign probe type.
    if rcp.split("P")[0] == "X":
        xp.append(True)
    else:
        xp.append(False)

    print("{}: {}  {:8.3f} - {:8.3f}  fg = {:.2f}".format(rcp, shot, time_start, time_end, avg_fg))

avg_machs = np.array(avg_machs)
fgs = np.array(list(fg.values()))
unf_xp = np.logical_and(unf, xp)
fav_xp = np.logical_and(~np.array(unf), xp)
unf_mp = np.logical_and(unf, ~np.array(xp))
fav_mp = np.logical_and(~np.array(unf), ~np.array(xp))

fig, ax = plt.subplots()
ax.scatter(fgs[unf_xp], avg_machs[unf_xp], color="tab:purple", marker="x", zorder=10, label="Unf. X-Point")
ax.scatter(fgs[fav_xp], avg_machs[fav_xp], color="tab:red", marker="x", zorder=10, label="Fav. X-Point")
ax.scatter(fgs[unf_mp], avg_machs[unf_mp], color="tab:purple", marker="^", zorder=10, label="Unf. MiMES")
ax.scatter(fgs[fav_mp], avg_machs[fav_mp], color="tab:red", marker="^", zorder=10, label="Fav. MiMES")
ax.axhline(0.0, linestyle="--", color="k")
ax.grid(zorder=1)
ax.legend(fontsize=10, ncol=2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlabel("Greenwald Fraction", fontsize=14)
ax.set_ylabel("Average Mach Number", fontsize=14)
fig.tight_layout()
fig.show()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharex=True, sharey=True)
ax1.scatter(fgs[unf_xp], avg_machs[unf_xp], color="tab:purple", marker="x", zorder=10, label="Unfavorable")
ax1.scatter(fgs[fav_xp], avg_machs[fav_xp], color="tab:red", marker="x", zorder=10, label="Favorable")
ax2.scatter(fgs[unf_mp], avg_machs[unf_mp], color="tab:purple", marker="^", zorder=10, label="Unfavorable")
ax2.scatter(fgs[fav_mp], avg_machs[fav_mp], color="tab:red", marker="^", zorder=10, label="Favorable")
ax1.set_title("XRCP", fontsize=16)
ax2.set_title("MRCP", fontsize=16)
ax1.axhline(0.0, linestyle="--", color="k")
ax1.grid(zorder=1)
ax1.legend(fontsize=12, ncol=1)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.set_xlabel("Greenwald Fraction", fontsize=14)
ax1.set_ylabel("Average Mach Number", fontsize=14)
ax2.axhline(0.0, linestyle="--", color="k")
ax2.grid(zorder=1)
ax2.legend(fontsize=12, ncol=1)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.set_xlabel("Greenwald Fraction", fontsize=14)
#ax2.set_ylabel("Average Mach Number", fontsize=14)
ax1.set_ylim([-0.5, 1.0])
fig.tight_layout()
fig.show()


df1 = dfs["MP184527_1"]
df2 = dfs["MP184527_2"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharex=True)
ax1.plot(df1["Rho"], df1["Te(eV)"], label="Plunge #1", color="r", lw=2)
ax1.plot(df2["Rho"], df2["Te(eV)"], label="Plunge #2", color="b", lw=2)
ax2.plot(df1["Rho"], df1["Ne (E18 m-3)"]*1e18, label="Plunge #1", color="r", lw=2)
ax2.plot(df2["Rho"], df2["Ne (E18 m-3)"]*1e18, label="Plunge #2", color="b", lw=2)
ax1.legend(fontsize=12)
ax2.legend(fontsize=12)
ax1.set_xlabel(r"$\rho$", fontsize=14)
ax2.set_xlabel(r"$\rho$", fontsize=14)
ax1.set_ylabel("Te (eV)", fontsize=14)
ax2.set_ylabel("ne (m-3)", fontsize=14)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax1.grid()
ax2.grid()
fig.tight_layout()
fig.show()
