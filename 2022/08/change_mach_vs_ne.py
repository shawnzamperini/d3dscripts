# Plot how the density changes vs a change in the Mach number (or flow).
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


path1 = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190485_1.tab"
path2 = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190485_2.tab"
df1 = pd.read_csv(path1, delimiter="\t")
df2 = pd.read_csv(path2, delimiter="\t")

fne1 = interp1d(df1["R(cm)"], df1["Ne(E18 m-3)"])
fne2 = interp1d(df2["R(cm)"], df2["Ne(E18 m-3)"])
#flow = "Vflow (km/s)"
flow = "Machn"
fm1 = interp1d(df1["R(cm)"], df1[flow])
fm2 = interp1d(df2["R(cm)"], df2[flow])

#comr = np.linspace(228, 233, 100)
comr = df1["R(cm)"][np.logical_and(df1["R(cm)"]>228, df1["R(cm)"]<233)]

fig, ax1 = plt.subplots(figsize=(5,4))

ax1.axhline(0, color="k")
ax11 = ax1.twinx()
ax1.plot(comr, fne2(comr)*1e18-fne1(comr)*1e18, lw=3, color="k")
ax11.plot(comr, fm2(comr)-fm1(comr), color="tab:green", lw=3)
ax1.set_xlabel("R (cm)", fontsize=14)
ax1.set_ylabel("Change in ne", fontsize=14)
ax11.set_ylabel("Change in Mach number", fontsize=14, color="tab:green")
ax11.tick_params(axis="y", labelcolor="tab:green")

fig.tight_layout()
fig.show()


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5,5), sharex=True)

ax1.axvline(228.9, color="k", linestyle="--")
ax2.axvline(228.9, color="k", linestyle="--")
ax2.axhline(0, color="k")
ax1.plot(df1["R(cm)"], df1["Ne(E18 m-3)"]*1e18, color="tab:red", lw=2)
ax1.plot(df2["R(cm)"], df2["Ne(E18 m-3)"]*1e18, color="tab:purple", lw=2)
ax2.plot(df1["R(cm)"], df1["Machn"], color="tab:red", lw=2)
ax2.plot(df2["R(cm)"], df2["Machn"], color="tab:purple", lw=2)

ax2.set_xlabel("R (cm)", fontsize=14)
ax1.set_ylabel("ne (m-3)", fontsize=14)
ax2.set_ylabel("Mach number", fontsize=14)
ax1.set_yscale("log")
ax1.grid(which="both", alpha=0.5)
ax2.grid(alpha=0.5)
ax1.set_ylim([5e17, 1e19])
ax1.set_xlim([226, 234])

fig.tight_layout()
fig.show()
