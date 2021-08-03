import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Some constants.
#ring = 40
mz_amu = 183.84
md_amu = 2.01
mz = mz_amu * 1.66e-27  # Mass of W in kg
md = md_amu * 1.66e-27   # Mass of D in kg
col_log = 15
charge = 8

#unf_path    = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-031a.nc"
fav_fl_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
fav_no_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
#unf = oedge_plots.OedgePlots(unf_path)
fav_no = oedge_plots.OedgePlots(fav_no_path)
fav_fl = oedge_plots.OedgePlots(fav_fl_path)

ring = 30
charge = 8
#s, ff = fav_fl.along_ring(ring, "ff", charge=charge, plot_it=False)
s, figf = fav_fl.along_ring(ring, "fig", charge=charge, plot_it=False)
s, fegf = fav_fl.along_ring(ring, "feg", charge=charge, plot_it=False)
s, fe = fav_fl.along_ring(ring, "fe", charge=charge, plot_it=False)
s, vz = fav_fl.along_ring(ring, "VELAVG", charge=charge, plot_it=False)
s, nz = fav_fl.along_ring(ring, "DDLIMS", charge=charge, plot_it=False)
s, ti = fav_fl.along_ring(ring, "KTIBS", plot_it=False)
s, ne = fav_fl.along_ring(ring, "KNBS", plot_it=False)
s, vi = fav_fl.along_ring(ring, "Velocity", plot_it=False)

# Can calculate FF with our returned vz.
tau_s = 1.47E13 * mz_amu * ti * np.sqrt(ti / md_amu) / \
        ((1 + md_amu / mz_amu) * ne * np.power(charge, 2) * col_log)
ff = mz * (vi - vz) / tau_s

gamma = vz * nz
fi = -mz / nz * np.gradient(vz * gamma, s)
fi = savgol_filter(fi, 21, 2)

fnet = ff + figf + fegf + fe + fi

fig, ax = plt.subplots()

ax.plot(s, ff, label="FF")
ax.plot(s, figf, label="FIG")
ax.plot(s, fegf, label="FEG")
ax.plot(s, fe, label="FE")
ax.plot(s, fi, label="FI")
ax.plot(s, fnet, color="k", label="Fnet")
ax.set_xlabel("S (m)")
ax.set_ylabel("Force (N)")
ax.legend()

fig.tight_layout()
fig.show()
