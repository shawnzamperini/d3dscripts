import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/Abrams_Plots_for_Jerome_062016.xlsx"
xl = pd.read_excel(path, sheet_name="167196")

x = xl["R-ROSP"].values
y = xl["Average"].values
yerr = np.abs(xl["20%*Avg"].values)

fig, ax1 = plt.subplots(figsize=(4, 3.5))
ax1.axvspan(142.5-142.68, 145.5-142.68, color="tab:cyan", alpha=0.25)
ax1.errorbar(x, y, yerr=yerr, lw=2, color="tab:cyan", elinewidth=1)
ax1.set_ylim([0, None])
ax1.set_xlim([None, 4])
ax1.set_yticks([0, 1e15, 2e15, 3e15, 4e15])
ax1.set_xlabel("Distance from Strike Point (cm)", fontsize=12)
ax1.set_ylabel(r"W Gross Erosion Rate ($\mathdefault{cm^{-2}s^{-1}}$)", fontsize=12)
ax1.text(0.78, 0.92, "#167196", transform=ax1.transAxes, fontsize=10)
ax1.text(0.3, 0.1, "W Ring", transform=ax1.transAxes, fontsize=14, color="tab:cyan")
fig.tight_layout()
fig.show()