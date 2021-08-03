import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Constants and values from Boedo 2001
rmrs = np.array([0.5, 5, 10])  # cm
vblob = np.array([2660, 1000, 330])  #m/s
tblob = 0.000015  # s
fblob = 3000
col_log = 15
charge = 15
eps0 = 8.85 * 10**(-12)  # F/m
elec = 1.609 * 10**(-19)

# Load the example RCP data.
xl_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2019/07/lp_with_a2.xlsx"
df = pd.read_excel(xl_path, sheet_name="Sheet3")

# Estimates of plasma values at each location.
rmrs_rcp = df["R-Rsep OMP Shifted (cm)"].values
te = df["Te (eV)"].values
ne = df["ne (1e18 m-3)"].values * 1e18
mD = 2.0 * 931.49e6 / (3e8)**2  # eV s2 / m2
mW = 74.0 * 931.49e6 / (3e8)**2
mD_amu = 2.0
mW_amu = 74.0
mW_kg = 74 * 1.66*10**(-27)
mD_kg = 2.0 * 1.66*10**(-27)
cs = np.sqrt(2 * te / mD)

# The stopping time.
tau_s = 1.47E13 * mW_amu * te * np.sqrt(te / mD_amu) / \
        ((1 + mD_amu / mW_amu) * ne * np.power(charge, 2) * col_log)

# Linearly interpolate alpha from 0 to 1 from LCFS to Wall
#alpha = np.linspace(0, 1, len(cs))
#alpha = 0.1 * np.exp(0.23 * rmrs_rcp)

# Normalize alpha where it is zero at R-Rsep=0, and 1 at 10.
rmrs_range = np.where(rmrs_rcp < 10)[0]
alpha = 1 - (tau_s - tau_s[rmrs_range].min()) / (tau_s[rmrs_range].max() - tau_s[rmrs_range].min())
#alpha = 0

# Estimate steps in connection length.
conn = np.zeros(len(rmrs_rcp))
conn[np.where(rmrs_rcp < 5)[0]] = 50
conn[np.where(np.logical_and(rmrs_rcp>=5, rmrs_rcp < 9))[0]] = 10
conn[np.where(rmrs_rcp > 9)[0]] = 1

conn = np.full(len(conn), 10)

# Estimate for time spent in SOL.
t_zpar = conn / cs

# Calculation of the collision frequency (Friedberg, Eq. 9.48).
vz = cs
vTi = cs
nu_zi = (1 / (4*np.pi) * elec**4 * ne / (eps0**2 * mW_kg * mD_kg) * 15) * 1 / (vz**3 + 1.3*vTi**3)

# Get linear estimates of vblob on the RCP domain.
r1 = np.where(rmrs_rcp <= 5)[0]
r2 = np.where(rmrs_rcp > 5)[0]
vblob1 = (1000 - 2660) / (5 - 0.5) * (rmrs_rcp[r1] - 0.5) + 2660
vblob2 = (330 - 1000) / (10 - 5) * (rmrs_rcp[r2] - 5) + 1000
vblob = np.append(vblob1, vblob2)

# Just a simple exponential to represent the decreasing influence of collisionality.
drag_coef = 25
drag_mod = np.exp(rmrs_rcp / drag_coef) - 1

drexb = vblob * (1 + drag_mod) * tblob * t_zpar * fblob
dr_fric = alpha * vblob * t_zpar
drrad = drexb + dr_fric

dperp_exb = drexb**2 / (2 * t_zpar)
dperp_fric = dr_fric**2 / (2 * t_zpar)
dperp = drrad ** 2 / (2 * t_zpar)
#dperp = dperp_exb + dperp_fric

"""
# Plotting.
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
ax1.plot(rmrs_rcp[r1], vblob1, "r-")
ax1.plot(rmrs_rcp[r2], vblob2, "b-")
ax1.set_ylabel("Blob radial velocity (m/s)")
ax2.plot(rmrs_rcp, tau_s)
ax2.set_ylabel("Tau_s (s)")
ax3.plot(rmrs_rcp, t_zpar)
ax3.set_ylabel("Connection length (m)")
ax4.plot(rmrs_rcp, dperp_exb, label="exb")
ax4.plot(rmrs_rcp, dperp_fric, label="fric")
ax4.plot(rmrs_rcp, dperp, label="sum")
ax4.set_ylabel("Dperp Contributions")
ax4.set_xlim([0, 10])
#ax4.set_ylim([0, 20])
#ax4.set_yscale("log")
ax4.legend()
fig.tight_layout()
fig.show()
"""

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharex=True)
#ax11 = ax1.twinx()
ax1.plot(rmrs_rcp[r1], vblob1, "k-", lw=4)
ax1.plot(rmrs_rcp[r2], vblob2, "k-", lw=4)
ax1.plot(rmrs_rcp[r1], vblob1*drag_mod[r1], "r-", lw=4)
ax1.plot(rmrs_rcp[r2], vblob2*drag_mod[r2], "r-", lw=4)
ax2.plot(rmrs_rcp, dperp_exb, "-", color="tab:red", lw=4)
#ax2.plot(10, 9, "*", color="tab:red", ms=20, mec="k", mew=2)
ax1.set_xlabel("R-Rsep (cm)", fontsize=16)
ax2.set_xlabel("R-Rsep (cm)", fontsize=16)
ax1.set_ylabel("Blob Radial Velocity (m/s)", fontsize=16)
ax2.set_ylabel("Corresponding Dperp Value", fontsize=16)
ax1.set_xlim([0, 11])
ax2.set_ylim([0, 10])
ax1.set_yticks(np.arange(0, 3000, 500))
ax1.set_yticklabels(np.arange(0, 3000, 500), fontsize=12)
ax2.set_yticks(np.arange(1, 11))
ax2.set_yticklabels(np.arange(1, 11), fontsize=12)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax1.grid()
ax2.grid()
ax1.set_ylim([0, 2500])
fig.tight_layout()
fig.show()
