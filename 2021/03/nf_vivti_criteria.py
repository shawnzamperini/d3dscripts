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
#fav_fl_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006d.nc"
#fav_no_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167277/d3d-167277-inj-006.nc"
#ring1 = 17
fav_no_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a.nc"
fav_fl_path = "/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/167247/d3d-167247-inj-034a.nc"
ring1 = 30

#unf = oedge_plots.OedgePlots(unf_path)
fav_no = oedge_plots.OedgePlots(fav_no_path)
fav_fl = oedge_plots.OedgePlots(fav_fl_path)

def calc_vinrt(vz, tau_s, nz, s):
    gamma = savgol_filter(vz, 21, 2) * savgol_filter(nz, 21, 2)
    vinrt = tau_s / nz * np.gradient(savgol_filter(vz * gamma, 21, 2), s)
    return vinrt

# Pick three regions of the plasma to compare the vz = vi + vTi approximation.
s1, vz_no1 = fav_no.along_ring(ring1, "VELAVG", charge=charge, plot_it=False)
s1, vz_fl1 = fav_fl.along_ring(ring1, "VELAVG", charge=charge, plot_it=False)
s1, nz_no1 = fav_no.along_ring(ring1, "DDLIMS", charge=charge, plot_it=False)
s1, nz_fl1 = fav_fl.along_ring(ring1, "DDLIMS", charge=charge, plot_it=False)
s1, vi_no1 = fav_no.along_ring(ring1, "Velocity", plot_it=False)
s1, vi_fl1 = fav_fl.along_ring(ring1, "Velocity", plot_it=False)
s1, ti1 = fav_no.along_ring(ring1, "KTIBS", plot_it=False)
s1, ne1 = fav_no.along_ring(ring1, "KNBS", plot_it=False)
s1, fig1 = fav_no.along_ring(ring1, "FIG", charge=charge, plot_it=False)
tau_s1 = 1.47E13 * mz_amu * ti1 * np.sqrt(ti1 / md_amu) / \
        ((1 + md_amu / mz_amu) * ne1 * np.power(charge, 2) * col_log)
vti1 = fig1 * tau_s1 / mz
vinrt_no1 = calc_vinrt(vz_no1, tau_s1, nz_no1, s1)
vinrt_fl1 = calc_vinrt(vz_fl1, tau_s1, nz_fl1, s1)

ring2 = 30
s2, vz_no2 = fav_no.along_ring(ring2, "VELAVG", charge=charge, plot_it=False)
s2, vz_fl2 = fav_fl.along_ring(ring2, "VELAVG", charge=charge, plot_it=False)
s2, nz_no2 = fav_no.along_ring(ring2, "DDLIMS", charge=charge, plot_it=False)
s2, nz_fl2 = fav_fl.along_ring(ring2, "DDLIMS", charge=charge, plot_it=False)
s2, vi_no2 = fav_no.along_ring(ring2, "Velocity", plot_it=False)
s2, vi_fl2 = fav_fl.along_ring(ring2, "Velocity", plot_it=False)
s2, ti2 = fav_no.along_ring(ring2, "KTIBS", plot_it=False)
s2, ne2 = fav_no.along_ring(ring2, "KNBS", plot_it=False)
s2, fig2 = fav_no.along_ring(ring2, "FIG", charge=charge, plot_it=False)
tau_s2 = 1.47E13 * mz_amu * ti2 * np.sqrt(ti2 / md_amu) / \
        ((1 + md_amu / mz_amu) * ne2 * np.power(charge, 2) * col_log)
vti2 = fig2 * tau_s2 / mz
vinrt_no2 = calc_vinrt(vz_no2, tau_s2, nz_no2, s2)
vinrt_fl2 = calc_vinrt(vz_fl2, tau_s2, nz_fl2, s2)

ring3 = 41
s3, vz_no3 = fav_no.along_ring(ring3, "VELAVG", charge=charge, plot_it=False)
s3, vz_fl3 = fav_fl.along_ring(ring3, "VELAVG", charge=charge, plot_it=False)
s3, nz_no3 = fav_no.along_ring(ring3, "DDLIMS", charge=charge, plot_it=False)
s3, nz_fl3 = fav_fl.along_ring(ring3, "DDLIMS", charge=charge, plot_it=False)
s3, vi_no3 = fav_no.along_ring(ring3, "Velocity", plot_it=False)
s3, vi_fl3 = fav_fl.along_ring(ring3, "Velocity", plot_it=False)
s3, ti3 = fav_no.along_ring(ring3, "KTIBS", plot_it=False)
s3, ne3 = fav_no.along_ring(ring3, "KNBS", plot_it=False)
s3, fig3 = fav_no.along_ring(ring3, "FIG", charge=charge, plot_it=False)
tau_s3 = 1.47E13 * mz_amu * ti3 * np.sqrt(ti3 / md_amu) / \
        ((1 + md_amu / mz_amu) * ne3 * np.power(charge, 2) * col_log)
vti3 = fig3 * tau_s3 / mz
vinrt_no3 = calc_vinrt(vz_no3, tau_s3, nz_no3, s3)
vinrt_fl3 = calc_vinrt(vz_fl3, tau_s3, nz_fl3, s3)


def plot_vel(ax, s, vz, vi, vti, vinrt, color):
    ax.plot(s, vi, linestyle="--", color=color)
    ax.plot(s, vti, linestyle=":", color=color)
    ax.plot(s, vinrt, linestyle="-.", color=color)
    ax.plot(s, vi+vti-vinrt, linestyle="-", color=color)
    ax.plot(s, vz, color="k")
    return ax

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))

ax1 = plot_vel(ax1, s1, vz_no1, vi_no1, vti1, vinrt_no1, color="r")
#ax1 = plot_vel(ax1, s1, vz_fl1, vi_fl1, vti1, vinrt_fl1, color="b")

ax2 = plot_vel(ax2, s2, vz_no2, vi_no2, vti2, vinrt_no2, color="r")
#ax2 = plot_vel(ax2, s2, vz_fl2, vi_fl2, vti2, vinrt_fl2, color="b")

ax3 = plot_vel(ax3, s3, vz_no3, vi_no3, vti3, vinrt_no3, color="r")
#ax3 = plot_vel(ax3, s3, vz_fl3, vi_fl3, vti3, vinrt_fl3, color="b")

fig.tight_layout()
fig.show()
