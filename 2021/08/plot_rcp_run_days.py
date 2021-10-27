# This is a simple script to just plot the RCP for each individual run day.
# Kinda just a way to get a first glance at the data, not really anything too
# specific in mind here yet.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# All the RCP data.
xl_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/" + \
  "rcp_master_detailed.xlsx"
master_xl = pd.ExcelFile(xl_path)

# Connection length data for a single shot from each day.
xl_path = "/Users/zamperini/Google Drive/My Drive/Research/Documents/2021/" + \
  "08/rcp_conns.xlsx"
conn_xl = pd.ExcelFile(xl_path)

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

def get_conn_data(conn_df, favorable, mrcp):
    if favorable:
        if mrcp:
            itf_x = conn_df["Psin"].dropna().values
            itf_y = conn_df["Lconn (m)"].dropna().values
            otf_x = conn_df["Psin.1"].dropna().values
            otf_y = conn_df["Lconn (m).1"].dropna().values
        else:
            itf_x = conn_df["Psin.2"].dropna().values
            itf_y = conn_df["Lconn (m).2"].dropna().values
            otf_x = conn_df["Psin.3"].dropna().values
            otf_y = conn_df["Lconn (m).3"].dropna().values
    else:
        if mrcp:
            itf_x = conn_df["Psin.1"].dropna().values
            itf_y = conn_df["Lconn (m).1"].dropna().values
            otf_x = conn_df["Psin"].dropna().values
            otf_y = conn_df["Lconn (m)"].dropna().values
        else:
            itf_x = conn_df["Psin.3"].dropna().values
            itf_y = conn_df["Lconn (m).3"].dropna().values
            otf_x = conn_df["Psin.2"].dropna().values
            otf_y = conn_df["Lconn (m).2"].dropna().values
    return {"itf_x":itf_x, "itf_y":itf_y, "otf_x":otf_x, "otf_y":otf_y}

def plot_conns(ax, conn_data):
    ax2 = ax.twinx()
    ax2.plot(conn_data["itf_x"], conn_data["itf_y"], color="k")
    ax2.plot(conn_data["otf_x"], conn_data["otf_y"], color="k", linestyle="--")

fig, axs = plt.subplots(5, 2, sharex=True, figsize=(10,7))
axs = axs.flatten()
conn_day_plotted = [False, False, False, False, False]
interps = {}
for sheet_name in master_xl.sheet_names:
    shot = int(sheet_name[2:8])
    ptype = sheet_name[:2]

    rcp_df = master_xl.parse(sheet_name)
    x = rcp_df["Psin"].values
    y = rcp_df["Mach"].values
    #y = rcp_df["Vflow (m/s)"].values
    #y = rcp_df["ne (1e18 m-3)"].values
    #y = rcp_df["Te (eV)"].values
    #if ptype == "MP":
    #    y = (rcp_df["Vf1 (V)"].values + rcp_df["Vf2 (V)"].values) / 2
    #else:
    #    y = rcp_df["Vf"].values

    if shot in day1:  # Unfavorable
        plot_num = 0
        if not conn_day_plotted[0]:
            conn_df = conn_xl.parse("184154", skiprows=1)
            xrcp_conn_data = get_conn_data(conn_df, False, False)
            mrcp_conn_data = get_conn_data(conn_df, False, True)
            plot_conns(axs[0], xrcp_conn_data)
            plot_conns(axs[1], mrcp_conn_data)
            f_xrcp_itf = interp1d(xrcp_conn_data["itf_x"], xrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_xrcp_otf = interp1d(xrcp_conn_data["otf_x"], xrcp_conn_data["otf_y"], fill_value="extrapolate")
            f_mrcp_itf = interp1d(mrcp_conn_data["itf_x"], mrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_mrcp_otf = interp1d(mrcp_conn_data["otf_x"], mrcp_conn_data["otf_y"], fill_value="extrapolate")
            interps["day1"] = {"f_xrcp_itf":f_xrcp_itf, "f_xrcp_otf":f_xrcp_otf,
              "f_mrcp_itf":f_mrcp_itf, "f_mrcp_otf":f_mrcp_otf}
            conn_day_plotted[0] = True

            if ptype == "XP":
                l_itf = f_xrcp_itf(x)
                l_otf = f_xrcp_otf(x)
                tot_conn = l_itf + l_otf
            else:
                l_itf = f_mrcp_itf(x)
                l_otf = f_mrcp_otf(x)
                tot_conn = l_itf + l_otf
                plot_num += 1
            l = l_itf - l_otf
            m = 2 / tot_conn - 1
            simp_m = m * l
            ax.plot(x, np.full(len(x), simp_m), color="r")

        f_day = interps["day1"]
    elif shot in day2:  # Favorable
        plot_num = 2
        conn_df = conn_xl.parse("184180", skiprows=1)
        if not conn_day_plotted[1]:
            conn_df = conn_xl.parse("184180", skiprows=1)
            xrcp_conn_data = get_conn_data(conn_df, True, False)
            mrcp_conn_data = get_conn_data(conn_df, True, True)
            plot_conns(axs[2], xrcp_conn_data)
            plot_conns(axs[3], mrcp_conn_data)
            f_xrcp_itf = interp1d(xrcp_conn_data["itf_x"], xrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_xrcp_otf = interp1d(xrcp_conn_data["otf_x"], xrcp_conn_data["otf_y"], fill_value="extrapolate")
            f_mrcp_itf = interp1d(mrcp_conn_data["itf_x"], mrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_mrcp_otf = interp1d(mrcp_conn_data["otf_x"], mrcp_conn_data["otf_y"], fill_value="extrapolate")
            interps["day2"] = {"f_xrcp_itf":f_xrcp_itf, "f_xrcp_otf":f_xrcp_otf,
              "f_mrcp_itf":f_mrcp_itf, "f_mrcp_otf":f_mrcp_otf}
            conn_day_plotted[1] = True

            if ptype == "XP":
                l_itf = f_xrcp_itf(x)
                l_otf = f_xrcp_otf(x)
                tot_conn = l_itf + l_otf
            else:
                l_itf = f_mrcp_itf(x)
                l_otf = f_mrcp_otf(x)
                tot_conn = l_itf + l_otf
                plot_num += 1
            l = l_itf - l_otf
            m = 2 / tot_conn - 1
            simp_m = m * l
            ax.plot(x, np.full(len(x), simp_m), color="r")

        f_day = interps["day2"]
    elif shot in day3:  # Favorable
        plot_num = 4
        conn_df = conn_xl.parse("184264", skiprows=1)
        if not conn_day_plotted[2]:
            conn_df = conn_xl.parse("184264", skiprows=1)
            xrcp_conn_data = get_conn_data(conn_df, True, False)
            mrcp_conn_data = get_conn_data(conn_df, True, True)
            plot_conns(axs[4], xrcp_conn_data)
            plot_conns(axs[5], mrcp_conn_data)
            f_xrcp_itf = interp1d(xrcp_conn_data["itf_x"], xrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_xrcp_otf = interp1d(xrcp_conn_data["otf_x"], xrcp_conn_data["otf_y"], fill_value="extrapolate")
            f_mrcp_itf = interp1d(mrcp_conn_data["itf_x"], mrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_mrcp_otf = interp1d(mrcp_conn_data["otf_x"], mrcp_conn_data["otf_y"], fill_value="extrapolate")
            interps["day3"] = {"f_xrcp_itf":f_xrcp_itf, "f_xrcp_otf":f_xrcp_otf,
              "f_mrcp_itf":f_mrcp_itf, "f_mrcp_otf":f_mrcp_otf}
            conn_day_plotted[2] = True

            if ptype == "XP":
                l_itf = f_xrcp_itf(x)
                l_otf = f_xrcp_otf(x)
                tot_conn = l_itf + l_otf
            else:
                l_itf = f_mrcp_itf(x)
                l_otf = f_mrcp_otf(x)
                tot_conn = l_itf + l_otf
                plot_num += 1
            l = l_itf - l_otf
            m = 2 / tot_conn - 1
            simp_m = m * l
            ax.plot(x, np.full(len(x), simp_m), color="r")

        f_day = interps["day3"]
    elif shot in day4:  # Unfavorable
        plot_num = 6
        conn_df = conn_xl.parse("184528", skiprows=1)
        if not conn_day_plotted[3]:
            conn_df = conn_xl.parse("184528", skiprows=1)
            xrcp_conn_data = get_conn_data(conn_df, False, False)
            mrcp_conn_data = get_conn_data(conn_df, False, True)
            plot_conns(axs[6], xrcp_conn_data)
            plot_conns(axs[7], mrcp_conn_data)
            f_xrcp_itf = interp1d(xrcp_conn_data["itf_x"], xrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_xrcp_otf = interp1d(xrcp_conn_data["otf_x"], xrcp_conn_data["otf_y"], fill_value="extrapolate")
            f_mrcp_itf = interp1d(mrcp_conn_data["itf_x"], mrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_mrcp_otf = interp1d(mrcp_conn_data["otf_x"], mrcp_conn_data["otf_y"], fill_value="extrapolate")
            interps["day4"] = {"f_xrcp_itf":f_xrcp_itf, "f_xrcp_otf":f_xrcp_otf,
              "f_mrcp_itf":f_mrcp_itf, "f_mrcp_otf":f_mrcp_otf}
            conn_day_plotted[3] = True

            if ptype == "XP":
                l_itf = f_xrcp_itf(x)
                l_otf = f_xrcp_otf(x)
                tot_conn = l_itf + l_otf
            else:
                l_itf = f_mrcp_itf(x)
                l_otf = f_mrcp_otf(x)
                tot_conn = l_itf + l_otf
                plot_num += 1
            l = l_itf - l_otf
            m = 2 / tot_conn - 1
            simp_m = m * l
            ax.plot(x, np.full(len(x), simp_m), color="r")

        f_day = interps["day4"]
    elif shot in day5:  # Unfavorable
        plot_num = 8
        conn_df = conn_xl.parse("187108", skiprows=1)
        if not conn_day_plotted[4]:
            conn_df = conn_xl.parse("187108", skiprows=1)
            xrcp_conn_data = get_conn_data(conn_df, False, False)
            mrcp_conn_data = get_conn_data(conn_df, False, True)
            plot_conns(axs[8], xrcp_conn_data)
            plot_conns(axs[9], mrcp_conn_data)
            f_xrcp_itf = interp1d(xrcp_conn_data["itf_x"], xrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_xrcp_otf = interp1d(xrcp_conn_data["otf_x"], xrcp_conn_data["otf_y"], fill_value="extrapolate")
            f_mrcp_itf = interp1d(mrcp_conn_data["itf_x"], mrcp_conn_data["itf_y"], fill_value="extrapolate")
            f_mrcp_otf = interp1d(mrcp_conn_data["otf_x"], mrcp_conn_data["otf_y"], fill_value="extrapolate")
            interps["day5"] = {"f_xrcp_itf":f_xrcp_itf, "f_xrcp_otf":f_xrcp_otf,
              "f_mrcp_itf":f_mrcp_itf, "f_mrcp_otf":f_mrcp_otf}
            conn_day_plotted[4] = True

            if ptype == "XP":
                l_itf = f_xrcp_itf(x)
                l_otf = f_xrcp_otf(x)
                tot_conn = l_itf + l_otf
            else:
                l_itf = f_mrcp_itf(x)
                l_otf = f_mrcp_otf(x)
                tot_conn = l_itf + l_otf
                plot_num += 1
            l = l_itf - l_otf
            m = 2 / tot_conn - 1
            simp_m = m * l
            ax.plot(x, np.full(len(x), simp_m), color="r")

        f_day = interps["day5"]

    # At this point the interpolation function should already have been
    # calculated.
    if ptype == "XP":
        l_itf = f_day["f_xrcp_itf"](x)
        l_otf = f_day["f_xrcp_otf"](x)
        tot_conn = l_itf + l_otf
    else:
        l_itf = f_day["f_mrcp_itf"](x)
        l_otf = f_day["f_mrcp_otf"](x)
        tot_conn = l_itf + l_otf
        plot_num += 1

    # This is a very crude assumption, but assume that at least near the
    # measurements the flows is ~L. So to get the true midpoint value of the
    # flows you would multiply by a factor proportional to how far from the
    # midpoint you are.
    l = l_itf - l_otf
    m = 2 / tot_conn - 1
    #delta_s = (tot_conn / 2) - l_itf
    #delta_mach = m * delta_s
    simp_m = m * l

    axs[plot_num].plot(x, y, color=shot_colors[shot], label=sheet_name)
    #axs[plot_num].plot(x, y+delta_mach, color=shot_colors[shot])

axs[0].set_xlim([1.0, 1.3])
for ax in axs:
    ax.set_ylim([-1.25, 1.25])
    #ax.set_ylim([-0.5, 0.5])
    #ax.set_ylim([0, 20])
    #ax.set_ylim([-60, 0])
    #ax.set_ylim([0, 60])
    #ax.set_ylim([0, 1000])
    ax.grid()
    ax.axhline(0, color="k")
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])

for ax in axs[[0, 2, 4, 6, 8]]:
    ax.legend(fontsize=5, ncol=3)

axs[0].text(0.8, 0.80, "Inner", transform=axs[0].transAxes, bbox=dict(fc="white"))
axs[0].text(0.8, 0.15, "Outer", transform=axs[0].transAxes, bbox=dict(fc="white"))

axs[0].set_title("X-Point RCP")
axs[1].set_title("Midplane RCP")

axs[0].set_ylabel("Unfavorable")
axs[2].set_ylabel("Favorable")
axs[4].set_ylabel("Favorable")
axs[6].set_ylabel("Unfavorable")
axs[8].set_ylabel("Unfavorable")

axs[8].set_xlabel("Psin")
axs[9].set_xlabel("Psin")
fig.tight_layout()
fig.show()
