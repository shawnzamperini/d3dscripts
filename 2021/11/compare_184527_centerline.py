import lim_plots
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import netCDF4


xlpath = "/Users/zamperini/My Drive/Research/Documents/2021/11/test_div_inj.xlsx"
itf = pd.read_excel(xlpath, sheet_name="MCP4 L1 Data")
otf = pd.read_excel(xlpath, sheet_name="MCP4 R1 Data")

machs = [0.0, 0.1, 0.2, 0.3, 0.4]
ncpaths = [
"/Users/zamperini/Documents/d3d_work/184527/mcp4-184527-tor240_m0_v3.nc",
"/Users/zamperini/Documents/d3d_work/184527/mcp4-184527-tor240_m01_v3.nc",
"/Users/zamperini/Documents/d3d_work/184527/mcp4-184527-tor240_m02_v3.nc",
"/Users/zamperini/Documents/d3d_work/184527/mcp4-184527-tor240_m03_v3.nc",
"/Users/zamperini/Documents/d3d_work/184527/mcp4-184527-tor240_m04_v3.nc"]
lps = []
for ncpath in ncpaths:
    lps.append(lim_plots.LimPlots(ncpath, combine_repeat_runs=False))

# Load in the simulation with a flat source as well.
ncflat = "/Users/zamperini/Documents/lim_runs/mcp4-184527-tor240_flatdif.nc"
lpflat = lim_plots.LimPlots(ncflat, combine_repeat_runs=False)

lams_itf_x     = itf["Axial Location [mm]"] / 10  # mm to cm
lams_itf_y     = itf["13C Excess"]
lams_itf_y_err = itf["13C Excess error"]
lams_otf_x     = otf["Axial Location (mm)"] / 10
lams_otf_y     = otf["13C Excess"]
lams_otf_y_err = otf["13C Excess error"]

def get_dep(ncpath):

    # Load netcdf file.
    nc = netCDF4.Dataset(ncpath)

    # Location of each P bin, and its width. Note syntax in pulling out data.
    ps     = nc.variables['PS'][:].data
    pwids  = nc.variables['PWIDS'][:].data

    # Array of poloidal locations (i.e. the center of each P bin).
    pol_locs = ps - pwids / 2.0

    # Distance cell centers along probe surface (i.e. the radial locations).
    rad_locs = nc.variables['ODOUTS'][:].data * 100  # m to cm

    # Get the centerline index (or closest to it).
    cline = np.abs(pol_locs).min()

    # This is the deposition array of the 2D collector probe faces.
    dep_arr = nc.variables['NERODS3'][0] * -1

    # Index the deposition array at the centerline for plotting.
    side1_x = rad_locs[np.where(rad_locs > 0.0)[0]]
    side1_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
    side2_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
    side2_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

    return {"side1_x":side1_x, "side1_y":side1_y, "side2_x":side2_x, "side2_y":side2_y}

deps = []
lim_itf_xs = []; lim_itf_ys = []
lim_otf_xs = []; lim_otf_ys = []
for i in range(0, len(ncpaths)):
    dep = get_dep(ncpaths[i])
    lim_itf_xs.append(dep["side2_x"])
    lim_itf_ys.append(dep["side2_y"])
    lim_otf_xs.append(dep["side1_x"])
    lim_otf_ys.append(dep["side1_y"])


# 184527 was Forward BT, which in the end means 3DLIM's ITF is actually the OTF
# side, and vice-versa.



#for i in range(0, len(lps)):
#    lim_dict = lps[i].centerline(show_plot=False)
#    lim_itf_xs.append(lim_dict["otf_x"] * 100)  # m to cm
#    lim_itf_ys.append(lim_dict["otf_y"])
#    lim_otf_xs.append(lim_dict["itf_x"] * 100)
#    lim_otf_ys.append(lim_dict["itf_y"])
flat_lim_dict = lpflat.centerline(show_plot=False)
flat_lim_itf_x = flat_lim_dict["otf_x"] * 100
flat_lim_itf_y = flat_lim_dict["otf_y"]
flat_lim_otf_x = flat_lim_dict["itf_x"] * 100
flat_lim_otf_y = flat_lim_dict["itf_y"]

window = 21
lams_itf_y_sg = savgol_filter(lams_itf_y, window, 3)[window:-window]
lams_otf_y_sg = savgol_filter(lams_otf_y, window, 3)[window:-window]
lams_itf_x_sg = lams_itf_x[window:-window]
lams_otf_x_sg = lams_otf_x[window:-window]

cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(0, 0.9, len(lps)))
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(lams_itf_x, lams_itf_y, alpha=0.3, color="tab:purple", lw=2)
ax1.plot(lams_otf_x, lams_otf_y, alpha=0.3, color="tab:red", lw=2)
ax1.plot(lams_itf_x_sg, lams_itf_y_sg, color="tab:purple", lw=3, label="ITF")
ax1.plot(lams_otf_x_sg, lams_otf_y_sg, color="tab:red", lw=3, label="OTF")
for i in range(0, 1):
    ax2.plot(lim_itf_xs[i], lim_itf_ys[i], color=colors[i], lw=3, label="M = {:.1f}".format(machs[i]))
    ax2.plot(lim_itf_xs[i], lim_itf_ys[i], color="tab:purple", lw=2)
    ax2.plot(lim_otf_xs[i], lim_otf_ys[i], color=colors[i], lw=3)
    ax2.plot(lim_otf_xs[i], lim_otf_ys[i], color="tab:red", lw=2)
ax1.set_ylim([-250, 3000])
ax2.set_ylim([-0.05, 0.5])
ax1.set_xlim([0, 6])
ax1.legend()
ax1.set_xlabel("Distance along probe (cm)", fontsize=14)
ax1.set_ylabel("LAMS Counts")
ax2.set_ylabel("3DLIM Deposition (arbitrary)")
fig.tight_layout()
fig.show()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
ax1.plot(lams_itf_x, lams_itf_y, alpha=0.3, color="tab:purple", lw=2)
ax1.plot(lams_itf_x_sg, lams_itf_y_sg, color="tab:purple", lw=3, label="ITF")
ax2.plot(lams_otf_x, lams_otf_y, alpha=0.3, color="tab:red", lw=2)
ax2.plot(lams_otf_x_sg, lams_otf_y_sg, color="tab:red", lw=3, label="OTF")
ax12 = ax1.twinx()
ax22 = ax2.twinx()
for i in range(0, len(lps)):
    ax12.plot(lim_itf_xs[i], lim_itf_ys[i], color="k", lw=3)
    ax12.plot(lim_itf_xs[i], lim_itf_ys[i], color=colors[i], lw=2, label="M = {:.1f}".format(machs[i]))
    ax22.plot(lim_otf_xs[i], lim_otf_ys[i], color="k", lw=3)
    ax22.plot(lim_otf_xs[i], lim_otf_ys[i], color=colors[i], lw=2, label="M = {:.1f}".format(machs[i]))
ax12.plot(flat_lim_itf_x, flat_lim_itf_y, color="k", lw=3)
ax12.plot(flat_lim_itf_x, flat_lim_itf_y, color="tab:cyan", lw=2, label="Flat")
ax22.plot(flat_lim_otf_x, flat_lim_otf_y, color="k", lw=3)
ax22.plot(flat_lim_otf_x, flat_lim_otf_y, color="tab:cyan", lw=2, label="Flat")
ax1.set_xlim([0, 6])
ax1.set_ylim([-250, 3000])
ax12.set_ylim([-0.05, 0.5])
ax22.set_ylim([-0.05, 0.5])
ax1.set_xlabel("Distance along probe (cm)", fontsize=14)
ax2.set_xlabel("Distance along probe (cm)", fontsize=14)
ax1.set_ylabel("LAMS Counts", fontsize=14)
ax22.set_ylabel("3DLIM Deposition (arbitrary)", fontsize=14)
ax22.legend()
ax1.set_title("ITF", fontsize=14)
ax2.set_title("OTF", fontsize=14)
fig.tight_layout()
fig.show()

fig, ax = plt.subplots()
for i in range(0, len(lps)):
    x = lps[i].nc.variables["divimp_probs"][:][0]
    y = lps[i].nc.variables["divimp_probs"][:][1]
    ax.plot(x, y, color=colors[i], lw=2, label="M = {:.1f}".format(machs[i]))
ax.legend()
ax.set_xlabel("Y (m)", fontsize=14)
ax.set_ylabel("Unnormalized Probability", fontsize=14)
fig.tight_layout()
fig.show()
