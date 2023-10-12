import oedge_plots
import pandas as pd
from scipy.interpolate import interp1d, griddata
import matplotlib.pyplot as plt
import pickle
import numpy as np


# Doesn't matter which, just need the grid and background solution.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-diff-002-predep.nc"
op = oedge_plots.OedgePlots(ncpath)

# Get RCP data for Mach numbers. Positive = Towards IT.
rcp = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP167195_2.tab", delimiter="\t")
# f_M = interp1d(rcp["R(cm)"] / 100, rcp["Machn"], bounds_error=False, fill_value=(rcp["Machn"].iloc[-1], rcp["Machn"].iloc[0]))

gfile_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/167195_4460.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]

rcp_shift = -0.015
rcp_coord = (rcp["R(cm)"].to_numpy() / 100 + rcp_shift, np.full(len(rcp["R(cm)"].values), -0.185))
rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), rcp_coord)

# Minus to align with DIVIMP direction.
f_M = interp1d(rcp_psin, -rcp["Machn"], bounds_error=False, fill_value=0)

# Get a profile along the RCP location. Need all the rings within this range.
div = op.fake_probe(2.20, 2.38, -0.188, -0.188, "Mach", rings_only=True)

# Go through one ring at a time and find what the corresponding Mach number was there.
new_m = []
for i in range(0, len(div["psin"])):
    new_m.append(f_M(div["psin"][i]))

fig, ax = plt.subplots(figsize=(5, 4))

ax.plot(div["psin"], div["Mach"], label="Old")
ax.plot(div["psin"], new_m, label="New")
ax.legend()

fig.tight_layout()
fig.show()

for i in range(0, len(new_m)):
    if div["ring"][i] >= 19 and div["ring"][i] <= 156:
        print("{}  {:.2f}".format(div["ring"][i], new_m[i]))