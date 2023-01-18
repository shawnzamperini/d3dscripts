# This script is a starting point. It takes one shot with RCP data and looks at the MAFOT field lines for each
# data point to see if the Mach number agrees with what it should be in regards to the simple SOL prescription (where
# Mach just linearly increases to the sound speed at each target.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve


# This script just focuses on looking into one shot. A later one can do multiple shots.
shot = 190440
rcp1 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190440_1.tab", delimiter="\t")
rcp2 = pd.read_csv("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190440_2.tab", delimiter="\t")

# Drop bad data up front.
rcp = rcp1.iloc[12:-3]
#rcp = rcp2.iloc[:-3]

# Load the MAFOT data.
#p1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp2000_conn_+1.dat"  # ITF
#m1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp2000_conn_-1.dat"  # OTF
p1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp3000_conn_+1.dat"  # ITF
m1_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp3000_conn_-1.dat"  # OTF
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
    "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
mafot1_p1 = pd.read_csv(p1_path, skiprows=52, names=columns, delimiter="\t")
mafot1_m1 = pd.read_csv(m1_path, skiprows=52, names=columns, delimiter="\t")

# Extract data for easy access.
rcp_r = rcp["R(cm)"].values
rcp_m = rcp["Machn"].values  # Positive = OTF - directed.
rcp_te = rcp["Te(eV)"].values
r1 = mafot1_p1["R (m)"] * 100
lconn1_itf = mafot1_p1["Lconn (km)"].values * 1000
lconn1_otf = mafot1_m1["Lconn (km)"].values * 1000

# Using just the connection length and the assumption that the Mach number linearly increases to 1 at each target,
# calculate what that M number would be based off the RCP data point location along the field lines.
f_itf = interp1d(r1, lconn1_itf)
f_otf = interp1d(r1, lconn1_otf)
rcp_lconn_itf = f_itf(rcp_r)
rcp_lconn_otf = f_otf(rcp_r)

# Treat the inner target as s=0.
rcp_s = rcp_lconn_itf
rcp_smax = rcp_lconn_otf + rcp_lconn_itf

# For each value, if we are less than smax/2, the Mach number is ramping towards the outer target (negative slope),
# and vice-versa if we are past smax/2.
slope = (1 + 1) / (rcp_smax - 0)
simple_m = slope * (rcp_s - 0) - 1

# More complex, we can solve the implicit equation for M derived for the SOL (1.42 in Peter's book). This requires
# knowing the coefficient C, which assuming diffusive cross-field transport, C = D/lambda_ne^2. We can find lamnda_ne
# right here.
rcp_ne = np.log(rcp["Ne(E18 m-3)"].values)

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(rcp_r, rcp_ne)

# Do a running window of lambda_ne at each location.
window_size = 5
lambda_ne = np.zeros(len(rcp_r))
for i in range(0, len(rcp_r)):
    window = np.full(len(rcp_r), False)
    if i < int(window_size / 2):
        window[:window_size] = True
    elif len(rcp_r) - i <= int(window_size / 2):
        window[-window_size:] = True
    else:
        window[i-int(window_size / 2):i+int(window_size / 2)+1] = True
    z = np.polyfit(rcp_r[window], rcp_ne[window], 1)
    p = np.poly1d(z)
    lambda_ne[i] = -1 / z[0]
    ax.plot(rcp_r[window], p(rcp_r[window]), color="k", lw=1)

ax.set_xlabel("R (cm)")
ax.set_ylabel("ln(ne)")
fig.tight_layout()
fig.show()

# For each RCP R value, solve the implicit equation for M.
cs = np.sqrt(2 * rcp_te / 931.49e6) * 3e8
C = (np.pi / 2 - 1) * cs / rcp_smax

def eq(M, s, C, cs):
    return 2 * np.arctan(M) - M - C * s / cs

D_solved = np.zeros(len(rcp_s))
vr_solved = np.zeros(len(rcp_s))
simple_m2 = np.zeros(len(rcp_r))
for i in range(0, len(rcp_r)):

    # This equation solves on a domain assuming [-L, L], and then solves for [0, L] assuming it is mirror negative on
    # for [-L, 0]. So we need to adjust the RCP s accordingly, and apply a negative if necessary.
    if rcp_s[i] > rcp_smax[i] / 2:
        s = rcp_s[i] - rcp_smax[i] / 2
        mult = 1.0
    else:
        s = rcp_smax[i] / 2 - rcp_s[i]
        mult = -1.0
    ans = fsolve(eq, 0.5, args=(s, C[i], cs[i]))[0] * mult
    simple_m2[i] = ans
    print("{:.2f} | {:.2f}: {:.2f}".format(s, rcp_smax[i], ans))

    # Using the same equation, let's extract a Dperp by using the measured M and backing out what C (and thus D) is.
    C_solved = (2 * np.arctan(rcp_m[i]) - rcp_m[i]) * cs[i] / s
    D_solved[i] = C_solved * np.square(lambda_ne[i]/100)

    # Similarly solve for a vr, here we assumed dv/dr ~ 0.
    vr_solved[i] = C_solved * lambda_ne[i] / 100

# What would be expected from lambda_ne.
D_solved2 = np.square(lambda_ne/100) * cs / (rcp_smax / 2)

# Load the blobby data.
blob1_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190440_1.tab"
blob2_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190440_2.tab"
blob = pd.read_csv(blob1_path, delimiter="\t")
blob_r = blob["R(cm)"]
avg_square_vr = blob["Vr(m/s)"] * blob["T_blob(e-6s)"] * 1e-6 * blob["Npeaks"] / 0.005  # 0.005 is the time period of the measurements.

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 7), sharex=True)

ax1.axhline(0, color="k", linestyle="--")
ax1.plot(rcp_r, simple_m, label="Linear")
ax1.plot(rcp_r, simple_m2, label="Diffusive")
ax1.plot(rcp_r, rcp_m, label="Measured")
ax2.plot(r1, lconn1_itf, label="ITF")
ax2.plot(r1, lconn1_otf, label="OTF")
ax3.plot(rcp_r, D_solved)
ax3.plot(rcp_r, D_solved2)
ax33 = ax3.twinx()
ax33.plot(rcp_r, vr_solved, color="r")
ax33.plot(blob_r, avg_square_vr, color="r", linestyle="--")

ax2.set_yscale("log")
ax1.set_xlim([224, 235])
ax2.set_ylim([3, 100])
ax1.legend()
ax2.legend()
ax1.set_ylabel("Mach")
ax3.set_xlabel("R (cm)")
ax2.set_ylabel("Connection Length (m)")
ax3.set_ylabel("Dperp (m2/s)")
ax33.set_ylabel("vr (m/s)", color="r")

fig.tight_layout()
fig.show()
