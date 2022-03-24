import pandas as pd
import matplotlib.pyplot as plt


path = "/Users/zamperini/My Drive/Research/Documents/2022/03/wall_coords.xlsx"
pol = pd.read_excel(path, sheet_name="Current", skiprows=1)
tor = pd.read_excel(path, sheet_name="Proposed", skiprows=1)

pol_r = pol["R (m)"].values
pol_z = pol["Z (m)"].values
tor_r = tor["R (m)"].values
tor_z = tor["Z (m)"].values

fig, ax = plt.subplots()
ax.plot(pol_r, pol_z, color="k")
ax.plot(tor_r, tor_z, color="k", linestyle="--")
ax.set_aspect("equal")
ax.set_axis_off()
fig.tight_layout()
fig.show()
