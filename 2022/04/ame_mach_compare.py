import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Century Gothic"
mimes = False

rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
if mimes:
    rcp267 = pd.read_excel(rcp_path, sheet_name="MP184267_1")
    rcp527 = pd.read_excel(rcp_path, sheet_name="MP184527_1")
    cut267 = 0.093886      #R = 2.354
    cut527 = 0.082037      #R = 2.344

else:
    rcp267 = pd.read_excel(rcp_path, sheet_name="XP184267_1")
    rcp527 = pd.read_excel(rcp_path, sheet_name="XP184527_1")[:-3]
    cut267 = 0.1913
    cut527 = 0.1913

xcol = "Dist. from sep."
c267 = "tab:red"
c527 = "tab:purple"


x267 = rcp267[xcol].values
y267 = rcp267["Mach"].values
x527 = rcp527[xcol].values
y527 = rcp527["Mach"].values

idx267 = np.where(x267 < cut267)[0]
idx527 = np.where(x527 < cut527)[0]

x267 = x267[idx267]
y267 = y267[idx267]
x527 = x527[idx527]
y527 = y527[idx527]

fig, ax = plt.subplots(figsize=(5, 4))

ax.plot(x267, y267, color=c267, lw=3)
ax.plot(x527, y527, color=c527, lw=3)
ax.axhline(0, color="k", linestyle="--")
ax.set_xlabel("Distance from separatrix (m)", fontsize=16)
ax.set_ylabel("Mach Number", fontsize=16)
ax.tick_params(which="both", labelsize=14)
ax.grid()

fig.tight_layout()
fig.show()
