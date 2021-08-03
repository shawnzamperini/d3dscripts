import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Choose what to plot.
#fav_shots = [184182, 184267, 184270]; fav_label = "Favorable"
#unf_shots = [184527, 184528, 184531]; unf_label = "Unfavorable"

#fav_shots = [184267]; fav_label = 184267
#unf_shots = [184527]; unf_label = 184527

fav_shots = [184267]; fav_label = "Favorable"
unf_shots = [184527]; unf_label = "Unfavorable"

#fav_shots = []; fav_label = "Favorable"
#unf_shots = [184527, 184528, 184531]; unf_label = "Unfavorable"

#fav_shots = [184182, 184267, 184270]; fav_label = "Favorable"
#unf_shots = []; unf_label = "Unfavorable"

# Either X-point of MiMES
ptype = "XP"
#ptype = "MP"

# Load data.
fav_dfs = []; unf_dfs = []
for shot in fav_shots:
    path1 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_1.tab".format(ptype, shot)
    path2 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_2.tab".format(ptype, shot)
    df1 = pd.read_csv(path1, sep="\t")
    df2 = pd.read_csv(path2, sep="\t")
    fav_dfs.append((df1, df2))
for shot in unf_shots:
    path1 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_1.tab".format(ptype, shot)
    path2 = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/{}{}_2.tab".format(ptype, shot)
    df1 = pd.read_csv(path1, sep="\t")
    df2 = pd.read_csv(path2, sep="\t")
    unf_dfs.append((df1, df2))
print("Files loaded.")

# Plotting.
fontsize = 16
lw = 3
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

if fav_shots == [] or unf_shots == []:
    fav_colors = ["tab:pink", "tab:cyan", "tab:brown"]
    unf_colors = ["tab:pink", "tab:cyan", "tab:brown"]
else:
    fav_colors = ["tab:red", "tab:red", "tab:red"]
    unf_colors = ["tab:purple", "tab:purple", "tab:purple"]

count = 0
for plunge1, plunge2 in fav_dfs:
    ax1.plot(plunge1["Rho"], plunge1["Ne (E18 m-3)"], color=fav_colors[count], lw=lw)
    ax2.plot(plunge1["Rho"], plunge1["Te(eV)"], color=fav_colors[count], lw=lw)
    ax3.plot(plunge1["Rho"], plunge1["Machn"], color=fav_colors[count], lw=lw)
    ax4.plot(plunge1["Rho"], plunge1["Vflow (km/s)"], color=fav_colors[count], lw=lw)
    ax1.plot(plunge2["Rho"], plunge2["Ne (E18 m-3)"], color=fav_colors[count], lw=lw)
    ax2.plot(plunge2["Rho"], plunge2["Te(eV)"], color=fav_colors[count], lw=lw)
    ax3.plot(plunge2["Rho"], plunge2["Machn"], color=fav_colors[count], lw=lw)
    ax4.plot(plunge2["Rho"], plunge2["Vflow (km/s)"], color=fav_colors[count], lw=lw)
    count += 1
count = 0
for plunge1, plunge2 in unf_dfs:
    ax1.plot(plunge1["Rho"], plunge1["Ne (E18 m-3)"], color=unf_colors[count], lw=lw)
    ax2.plot(plunge1["Rho"], plunge1["Te(eV)"], color=unf_colors[count], lw=lw)
    ax3.plot(plunge1["Rho"], plunge1["Machn"], color=unf_colors[count], lw=lw)
    ax4.plot(plunge1["Rho"], plunge1["Vflow (km/s)"], color=unf_colors[count], lw=lw)
    ax1.plot(plunge2["Rho"], plunge2["Ne (E18 m-3)"], color=unf_colors[count], lw=lw)
    ax2.plot(plunge2["Rho"], plunge2["Te(eV)"], color=unf_colors[count], lw=lw)
    ax3.plot(plunge2["Rho"], plunge2["Machn"], color=unf_colors[count], lw=lw)
    ax4.plot(plunge2["Rho"], plunge2["Vflow (km/s)"], color=unf_colors[count], lw=lw)
    count += 1
#ax1.set_ylim([0, 12])
#ax2.set_ylim([0, 40])
ax1.grid(which="both")
ax2.grid(which="both")
ax3.grid()
ax4.grid()
ax3.set_ylim([-1, 1])
ax4.set_ylim([-25, 25])
ax1.set_yscale("log")
ax2.set_yscale("log")
ax3.axhline(0.0, color="k", linestyle="--")
ax4.axhline(0.0, color="k", linestyle="--")
ax3.set_xlabel("Rho", fontsize=fontsize)
ax4.set_xlabel("Rho", fontsize=fontsize)
ax1.set_ylabel("ne (1e18 m-3)", fontsize=fontsize)
ax2.set_ylabel("Te (eV)", fontsize=fontsize)
ax3.set_ylabel("Mach", fontsize=fontsize)
ax4.set_ylabel("Velocity (km/s)", fontsize=fontsize)

if fav_shots == []:
    custom_lines = [Line2D([0], [0], color=unf_colors[0]),
                    Line2D([0], [0], color=unf_colors[1]),
                    Line2D([0], [0], color=unf_colors[2])]
    ax2.legend(custom_lines, unf_shots, fontsize=15)
    ax3.legend(custom_lines, unf_shots, fontsize=15)
elif unf_shots == []:
    custom_lines = [Line2D([0], [0], color=fav_colors[0]),
                    Line2D([0], [0], color=fav_colors[1]),
                    Line2D([0], [0], color=fav_colors[2])]
    ax2.legend(custom_lines, fav_shots, fontsize=15)
    ax3.legend(custom_lines, fav_shots, fontsize=15)
else:
    custom_lines = [Line2D([0], [0], color="tab:red"),
                    Line2D([0], [0], color="tab:purple")]
    ax2.legend(custom_lines, [fav_label, unf_label], fontsize=12)
    ax3.legend(custom_lines, [fav_label, unf_label], fontsize=12)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)
fig.tight_layout()
fig.show()
