# Point of this script is simply to interpolate connection length data onto
# a 3DLIM set of R cell coordinates.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


conn_path = "/Users/zamperini/My Drive/School/Tennessee/Research/"+\
  "Collector Probe Excel Sheets/Connection Lengths/167247/167247.xlsx"
conn_df = pd.read_excel(conn_path, sheet_name="MAFOT ITF", skiprows=2)

# Tip of CP for this shot in cm (A8).
ptip = 227.787

lim_r = -(conn_df["R (m)"].values - ptip / 100)
conns = conn_df["Connection Length (km)"].values * 1000
f_L = interp1d(lim_r, conns, bounds_error=False, fill_value=conns[-1])

fig, ax = plt.subplots()
ax.plot(lim_r, conn_df["Connection Length (km)"] * 1000)


input_r = np.arange(-0.1475, 0.0225, 0.0025)
input_l = f_L(input_r)

ax.plot(input_r, input_l)

ax.grid()
ax.set_xlabel("3DLIM R")
ax.set_ylabel("Connection Length (m)")
fig.tight_layout()
fig.show()

print("Connection lengths")
for i in input_l:
    print(i)
