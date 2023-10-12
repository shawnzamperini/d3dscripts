import get_lp
import MDSplus
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from gadata import gadata


# Estimation of the power into SOL during 167196. PINJ - PRAD_CORE
psol = 2.6e6
psol_min_prad_tot = 1.5e6

# Load data during the strike point sweep.
lps = get_lp.plot_lps(shot=167195, tmin=4000, tmax=5000, bins=25, tunnel=False)

# Load strike point location so we can create a pseudo-R of every point during sweep.
conn = MDSplus.Connection("atlas.gat.com")
rvsout = gadata("RVSOUT", 167195, connection=conn)
f_rvsout = interp1d(rvsout.xdata, rvsout.zdata)

# Extract data from the relevant probes.
keep_pnames = ["probe 23", "probe 25", "probe 31", "probe 33", "probe 35"]
pmask = [True if p in keep_pnames else False for p in lps["pnames"]]
r = np.array(lps["R"])[pmask]
rminrsep = np.array(lps["rminrsep"])[pmask]
qpar = np.array(lps["heatflux (W/cm2)"])[pmask] * 10000  # W/cm2 --> W/m2
times = np.array(lps["avg_time"])[pmask]

# Calculate what the "psuedo-R" is of each data point when mapping back to a stationary equilibrium. Sort values
# by this data.
pseudo_r = f_rvsout(times) + rminrsep
sort_idx = np.argsort(pseudo_r)
pseudo_r = pseudo_r[sort_idx]
r = r[sort_idx]
rminrsep = rminrsep[sort_idx]
qpar = qpar[sort_idx]
times = times[sort_idx]

# Integrate to get total power to target.
qpar_1d = qpar * 2 * np.pi * pseudo_r  # Multiply by toroidal length, thus W/m now.
run_int = [np.trapz(qpar_1d[:i], pseudo_r[:i]) for i in range(1, len(qpar_1d))]

# Plotting.
fig, ax1 = plt.subplots(figsize=(5, 4))
ax2 = ax1.twinx()
ax1.scatter(pseudo_r, qpar)
ax2.plot(pseudo_r[1:], run_int, color="tab:red")
ax2.axhline(psol, color="tab:red", linestyle="--")
ax2.axhline(psol_min_prad_tot, color="tab:red", linestyle="-.")
ax1.set_xlabel("R (m)")
ax1.set_ylabel("qpar (W/m2)")
ax2.set_ylabel("Integrated Power (W)", color="tab:red")
fig.tight_layout()
fig.show()
