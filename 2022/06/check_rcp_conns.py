# Script to plot the RCP data with the connection length to see if it can
# explain different Mach profiles for 167192 and 167195.
import pandas as pd
import matplotlib.pyplot as plt


# Load the RCP data.
xlpath = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_167192-195.xlsx"
dfs192 = []
for plunge in ["MP167192_1", "MP167192_2"]:
    dfs192.append(pd.read_excel(xlpath, sheet_name=plunge))
dfs195 = []
for plunge in ["MP167195_1", "MP167195_2"]:
    dfs195.append(pd.read_excel(xlpath, sheet_name=plunge))

# Load the MAFOT data.
mafot_dat = {}
root = "/Users/zamperini/Documents/d3d_work/mafot_files/"
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
for shot in [167192, 167195]:
    tmp_dict = {}
    for dir in ["p1", "m1", "both"]:
        path = "{}{}/lam_efitwall_{}.dat".format(root, shot, dir)
        df = pd.read_csv(path, skiprows=52, names=columns, delimiter="\t")
        tmp_dict[dir] = df
    mafot_dat[shot] = tmp_dict

# Plots.
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharex=True)

ax1.axvline(2.283, color="k", linestyle="--")
ax2.axvline(2.296, color="k", linestyle="--")
ax1.axhline(0, color="k", linestyle="--")
ax2.axhline(0, color="k", linestyle="--")
ax1.plot(dfs192[0]["R(cm)"]/100, dfs192[0]["Machn"], color="r", marker="o")
ax2.plot(dfs195[0]["R(cm)"]/100, dfs195[0]["Machn"], color="r", marker="o")
ax1.plot(dfs192[1]["R(cm)"]/100, dfs192[1]["Machn"], color="r", marker="o")
ax2.plot(dfs195[1]["R(cm)"]/100, dfs195[1]["Machn"], color="r", marker="o")
ax1.set_xlim([2.25, 2.37])
ax1.set_ylim([-1, 1])
ax2.set_ylim([-1, 1])

ax11 = ax1.twinx()
ax22 = ax2.twinx()
ax11.plot(mafot_dat[167192]["both"]["R (m)"], mafot_dat[167192]["both"]["Lconn (km)"]*1000, color="k")
ax22.plot(mafot_dat[167195]["both"]["R (m)"], mafot_dat[167195]["both"]["Lconn (km)"]*1000, color="k")
ax11.set_yscale("log")
ax22.set_yscale("log")
ax11.set_ylim(1, 100)
ax22.set_ylim(1, 100)

ax1.set_title(167192)
ax2.set_title(167195)
ax1.set_xlabel("R (m)")
ax2.set_xlabel("R (m)")
ax1.set_ylabel("Mach", color="r")
ax22.set_ylabel("Connection Length (m)")
ax1.tick_params(axis='y', labelcolor="r")
ax2.tick_params(axis='y', labelcolor="r")
ax1.text(2.287, 0.5, "Helicon")

fig.tight_layout()
fig.show()

fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(9, 4), sharex=True, sharey=True)
ax3.axvline(2.283, color="b", linestyle="--")
ax4.axvline(2.296, color="b", linestyle="--")
ax3.plot(mafot_dat[167192]["p1"]["R (m)"], mafot_dat[167192]["p1"]["Lconn (km)"]*1000, color="r", label="OTF")
ax4.plot(mafot_dat[167195]["p1"]["R (m)"], mafot_dat[167195]["p1"]["Lconn (km)"]*1000, color="r", label="OTF")
ax3.plot(mafot_dat[167192]["m1"]["R (m)"], mafot_dat[167192]["m1"]["Lconn (km)"]*1000, color="b", label="ITF")
ax4.plot(mafot_dat[167195]["m1"]["R (m)"], mafot_dat[167195]["m1"]["Lconn (km)"]*1000, color="b", label="ITF")
ax3.set_yscale("log")
ax4.set_yscale("log")
ax3.set_ylim(1, 100)
ax4.set_ylim(1, 100)
ax3.set_title(167192)
ax4.set_title(167195)
ax3.set_ylabel("Connection Length (m)")
ax3.set_xlabel("R (m)")
ax4.set_xlabel("R (m)")
ax3.text(2.287, 40, "Wall Structure")
ax3.legend()
ax3.set_xlim([2.25, 2.37])
fig.tight_layout()
fig.show()
