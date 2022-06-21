# Plot the plunges from RCP data from 167192-195
import matplotlib.pyplot as plt
import pandas as pd


xlpath = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_167192-195.xlsx"
plunges = ["MP167192_1", "MP167192_2", "MP167193_1", "MP167193_2", "MP167194_1",
    "MP167194_2", "MP167195_1", "MP167195_2"]
dfs = []
for plunge in plunges:
    dfs.append(pd.read_excel(xlpath, sheet_name=plunge))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 7))
ax3.axhline(0.0, color="k", linestyle="--")
for i in range(0, len(dfs)):
    df = dfs[i]
    r = df["R(cm)"] / 100
    te = df["Te(eV)"]
    ne = df["Ne(E18 m-3)"] * 1e18
    m = df["Machn"]
    isat = df["Isat(A)"]
    ax1.plot(r, te, label=plunges[i])
    ax2.plot(r, ne)
    ax3.plot(r, m)
    ax4.plot(r, isat)

fig.supxlabel("R (m)")
ax1.set_ylabel("Te (eV)")
ax2.set_ylabel("ne (m-3)")
ax3.set_ylabel("Mach number")
ax4.set_ylabel("Isat (A)")
ax1.legend()
ax1.set_ylim([0, 60])
ax2.set_ylim([0, 2e19])
ax3.set_ylim([-1, 1])

fig.tight_layout()
fig.show()
