# A 2D plot comparing the inertail volcity comparison.
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

fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
#fav_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a.nc"
fav = oedge_plots.OedgePlots(fav_path)


# Gather 2D plasma parameters and velocities.
ti = fav.read_data_2d("KTIBS")
ne = fav.read_data_2d("KNBS")
fig = fav.calculate_forces("FIG", charge=charge)
vz = fav.read_data_2d("VELAVG", charge=charge)
nz = fav.read_data_2d("DDLIMS", charge=charge)
vi = fav.read_data_2d_kvhs_t13()
s = fav.read_data_2d("KSS")

tau_s = 1.47E13 * mz_amu * ti * np.sqrt(ti / md_amu) / \
        ((1 + md_amu / mz_amu) * ne * np.power(1, 2) * col_log)
vti = fig * tau_s / mz
vinrt = tau_s / nz * np.gradient(vz * vz * nz, s)

comp1 = (vz - (vi + vz)) / vz
comp2 = (vz - (vi + vz - vinrt)) / vz
fav.plot_contour_polygon("KTEBS", own_data=comp1, vmin=-1.0, vmax=1.0, normtype="symlin")
fav.plot_contour_polygon("KTEBS", own_data=comp2, vmin=-1.0, vmax=1.0, normtype="symlin")
