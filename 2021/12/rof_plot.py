import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


a2_path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2_path)

a2_itf_x = a2["Distance from Tip D (cm)"].values
a2_otf_x = a2["Distance from Tip U (cm)"].values
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values

cmap = plt.get_cmap("magma")
colors = cmap(np.linspace(0, 0.9, 4))

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(a2_itf_x, a2_itf_y,     color=colors[0], lw=3)
ax.plot(a2_itf_x, a2_itf_y*1.5, color=colors[1], lw=3)
ax.plot(a2_itf_x, a2_itf_y*2.0, color=colors[2], lw=3)
ax.plot(a2_itf_x, a2_itf_y*2.5, color=colors[3], lw=3)
ax.set_xlabel("Distance along probe (cm)", fontsize=14)
ax.set_ylabel("W183 Areal Depostion (W/cm2)", fontsize=14)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.tick_params(which="both", labelsize=11)
ax.set_xlim([0, 6])
fig.tight_layout()
fig.show()

x = [0.23, 0.54, 0.77, 0.96, 1.22]
y = [0.27, 0.49, 0.75, 0.94, 1.24]
y2 = [0.21, 0.56, 0.82, 0.91, 1.14]

fig, ax = plt.subplots()

ax.scatter(x, y, color="tab:purple")

ax2 = ax.twinx()
ax2.scatter(x, y2, color="tab:red")

ax.set_xlabel("Injected W183", fontsize=14)
ax.set_ylabel("Equation", fontsize=14)
ax2.set_ylabel("SXR Measurement", fontsize=14)

fig.tight_layout()
fig.show()
