# Script to estimate a corresponding diffusion coefficient from some real
# measurements of Te and connections lengths.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


# Load RCP data.
xl_path = "/Users/zamperini/My Drive/Research/Documents/2021/11/lp_167196.xlsx"
rcp = pd.read_excel(xl_path, sheet_name="Data")
r = rcp["R (cm)"].values / 100  # cm to m
te = rcp["Te (eV)"].values
ne = rcp["ne (1e18 m-3)"].values

# Bin into 1 cm bins, average.
bin_mins = np.arange(r.min(), r.max(), 0.01)
bin_idxs = np.digitize(r, bin_mins)
r_binned  = []
te_binned = []
ne_binned = []
for b in np.unique(bin_idxs):
    idx = bin_idxs == b
    r_binned.append(np.nanmean(r[idx]))
    te_binned.append(np.nanmean(te[idx]))
    ne_binned.append(np.nanmean(ne[idx]))
r_binned = np.array(r_binned)
te_binned = np.array(te_binned)
ne_binned = np.array(ne_binned)

# Replace with interpolated functions to get a simpler R base.
ft = interp1d(r_binned, te_binned, bounds_error=False, fill_value=(te_binned[0], te_binned[-1]))
fn = interp1d(r_binned, ne_binned, bounds_error=False, fill_value=(ne_binned[0], ne_binned[-1]))
r_binned = np.arange(r.min(), r.max(), 0.01)
te_binned = ft(r_binned)
ne_binned = fn(r_binned)

# Load connection lengths in each direction.
mafot_file1 = "/Users/zamperini/Documents/d3d_work/167196/lam_conns_-1.dat"
mafot_file2 = "/Users/zamperini/Documents/d3d_work/167196/lam_conns_+1.dat"
print("Loading MAFOT runs...")
columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
  "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
try:
    conn1 = pd.read_csv(mafot_file1, skiprows=52, names=columns,
      delimiter="\t")
except FileNotFoundError:
    print("Error: Unable to find file: {}".format(mafot_file1))
    print("Exiting")
    sys.exit()
try:
    conn2 = pd.read_csv(mafot_file2, skiprows=52, names=columns,
      delimiter="\t")
except FileNotFoundError:
    print("Error: Unable to find file: {}".format(mafot_file2))
    print("Exiting")
    sys.exit()

conn_r = conn1["R (m)"].values
conn_tot = (conn1["Lconn (km)"].values + conn1["Lconn (km)"].values) * 1000  # km to m

# At each RCP bin, estmate the characteristic parallel time.
mD = 2.0 * 931.49e6 / (3e8)**2  # eV s2 / m2
cs = np.sqrt(2 * te_binned / mD)
tzpar = []
nearest_conns = []
for i in range(0, len(cs)):
    near_conn_idx = np.where(np.abs(r_binned[i] - conn_r) == np.abs(r_binned[i] - conn_r).min())[0][0]
    tmp_tot = conn_tot[near_conn_idx]
    tzpar.append((tmp_tot/2) / cs[i])
    nearest_conns.append(tmp_tot)
tzpar = np.array(tzpar)
nearest_conns = np.array(nearest_conns)

# Interpolation functions to get some real rough estimates of blob characteristics.
rsep = 2.217
r_blob = rsep + np.array([0.005, 0.05, 0.10])
e_blob = np.array([4000, 1500, 500])
t_blob = 1.5e-5
v_blob = np.array([2660, 1000, 330])
fv_blob = interp1d(r_blob, v_blob, bounds_error=False, fill_value=(v_blob.max(), v_blob.min()))

# A real rough hypothetical Dperp.
#dperp = np.array([0.3, 1, 10])
#f_dperp = interp1d(r_blob, dperp, bounds_error=False, fill_value=(dperp.min(), dperp.max()))

# Between each binned R value (deltaR), see what
dr = np.array([r_binned[i] - r_binned[i-1] for i in range(1, len(r_binned))])
dt = dr / fv_blob(r_binned[1:])
#dperp = np.square(dr) / (2 * tzpar[1:])
dperp = np.square(dr) / (2 * dt)


#fig, (ax, ax3) = plt.subplots(2, 1, sharex=True)
fig, ax = plt.subplots(figsize=(5,4))
ax2 = ax.twinx()
ax2.plot(r_binned[1:], savgol_filter(dperp, 3, 1), color="tab:red", lw=2)
ax.plot(r_binned, savgol_filter(fv_blob(r_binned), 3, 1), color="tab:purple", lw=2)
ax2.set_title("#167196", fontsize=18)
ax.set_xlabel("R (m)", fontsize=14)
ax2.set_ylabel(r"$\mathdefault{D_{rad}}$", fontsize=14, color="tab:red")
ax.set_ylabel("Blob Velocity (m/s)", fontsize=14, color="tab:purple")
ax2.spines["top"].set_visible(False)
ax2.spines["left"].set_color("tab:purple")
ax2.spines["right"].set_color("tab:red")
ax.tick_params(axis="y", colors="tab:purple")
ax2.tick_params(axis="y", colors="tab:red")
fig.tight_layout()
fig.show()

fig, ax3 = plt.subplots()
ax4 = ax3.twinx()
ax3.plot(conn_r, conn_tot, color="tab:cyan", lw=2)
ax4.plot(r_binned, ne_binned, color="tab:pink", lw=2, zorder=40)
ax4.scatter(r_binned, ne_binned, marker="^", color="tab:pink", edgecolors="k", zorder=50)
ax4.set_yscale("log")
ax3.set_yscale("log")
ax3.set_xlabel("R (m)", fontsize=14)
ax3.set_ylabel("Lconn (m)", fontsize=14, color="tab:cyan")
ax4.set_ylabel("ne (1e18 m-3)", fontsize=14, color="tab:pink")
ax3.set_xlim(2.24, 2.38)
ax3.set_ylim([0.1, 200])
ax4.spines["left"].set_color("tab:cyan")
ax4.spines["right"].set_color("tab:pink")
ax3.tick_params(axis="y", which="both", colors="tab:cyan")
ax4.tick_params(axis="y", which="both", colors="tab:pink")

fig.tight_layout()
fig.show()
