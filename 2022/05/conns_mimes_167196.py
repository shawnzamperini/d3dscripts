import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit


# Load MAFOT data.
shot = 184267
path1 = "/Users/zamperini/Documents/d3d_work/mafot_files/{}/lam_mimes_plunge.dat".format(shot)
path2 = "/Users/zamperini/Documents/d3d_work/mafot_files/{}/lam_mimes_plunge_tor.dat".format(shot)
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
df1 = pd.read_csv(path1, skiprows=52, names=columns, delimiter="\t")
df2 = pd.read_csv(path2, skiprows=52, names=columns, delimiter="\t")

# Load RCP data.
if shot == 167196:
    path3 = "/Users/zamperini/My Drive/Research/Documents/2022/05/conn_comparison_167196.xlsx"
    df3 = pd.read_excel(path3, sheet_name="For Import")[["R (cm)", "R-Rsep OMP (m)", "Te (eV)", "ne (1e18 m-3)"]].dropna()
    rcp_r = df3["R (cm)"].values
    rcp_ne = df3["ne (1e18 m-3)"].values
    split_r = 233
    width = 0.75
elif shot == 184267:
    path3 = "/Users/zamperini/My Drive/Research/Data/rcp_data/rcp_master_detailed.xlsx"
    df3 = pd.read_excel(path3, sheet_name="MP184267_1")[["R (cm)","ne (1e18 m-3)"]].dropna()
    rcp_r = df3["R (cm)"].values
    rcp_ne = df3["ne (1e18 m-3)"].values
    split_r = 232
    width = 1

# Do exponential fits in each region.
def exp_fit(x, a, b):
    return a * np.exp(-x * b)

region1 = rcp_r < split_r
region2 = rcp_r > split_r

# More reasonable X values so fitting goes smoothly.
x1 = rcp_r[region1] - rcp_r[region1].min()
x2 = rcp_r[region2] - rcp_r[region2].min()
n1 = rcp_ne[region1]
n2 = rcp_ne[region2]
popt1, pcov1 = curve_fit(exp_fit, x1, n1, maxfev=5000, p0=(10, 1))
popt2, pcov2 = curve_fit(exp_fit, x2, n2, maxfev=5000, p0=(10, 10))

# Convert back to normal R coordinates, calculate fitted values.
r1 = x1 + rcp_r[region1].min()
r2 = x2 + rcp_r[region2].min()
n1f = exp_fit(x1, *popt1)
n2f = exp_fit(x2, *popt2)

print("lambda1 = {:.3f}".format(1/popt1[1]))
print("lambda2 = {:.3f}".format(1/popt2[1]))
print("lambda2/lambda1 = {:.2f}".format(popt1[1] / popt2[1]))

# I know there are 500 R points at one Z location.
r = df1["R (m)"].unique() * 100
z = np.full(len(r), -0.188)
l1 = df1["Lconn (km)"].values
l2 = df2["Lconn (km)"].values

# Use RCP R-Rsep OMP to R mapping to convert the MAFOT R values to OMP coords.
#normal_r = df3["R (cm)"].values
#z = np.polyfit(normal_r, rcp_r, 1)
#p = np.poly1d(z)
#r_mf = p(r*100)   # R (m) to R-Rsep OMP (m)

# Estimate change in Lconn by going half a cm back anf then forward.
split_idx = np.where(np.abs(r-split_r) == np.abs(r-split_r).min())[0][0]
left = np.where(np.abs(r-r[split_idx]-width) == np.abs(r-r[split_idx]-width).min())[0][0]
right = np.where(np.abs(r-r[split_idx]+width) == np.abs(r-r[split_idx]+width).min())[0][0]
print("Lconn change 1: {:.3f} m".format((l1[left] - l1[right]) * 1000))
print("Lconn change 2: {:.3f} m".format((l2[left] - l2[right]) * 1000))

# Plotting.
plt.rcParams["font.family"] = "Century Gothic"
fig, ax1 = plt.subplots()
ax11 = ax1.twinx()

ax1.plot(r, l1, color="r", lw=3)
ax1.plot(r, l2, linestyle="--", color="r", lw=3)
ax1.axhline(l1[left], color="r")
ax1.axhline(l1[right], color="r")

ax11.scatter(rcp_r, rcp_ne, c="k", alpha=0.4, linewidths=0, s=15)
ax11.plot(r1, n1f, color="k", lw=3)
ax11.plot(r2, n2f, color="k", lw=3)

ax1.axvline(split_r, color="k", linestyle="--")

ax1.set_ylim([0.0005, 0.1])
ax1.set_xlim([222, 240])
ax1.set_xlabel("R (cm)")
ax1.set_ylabel("Connection Length (km)", color="r")
ax11.set_ylabel("ne (1e18 m-3)")
ax11.set_yscale("log")
ax1.set_yscale("log")
fig.tight_layout()
fig.show()
