# Script to compare a flat source with the DIVIMP source to A2 deposition,
# indicating that connection length alone cannot account for the deposition
# patterns (ITF > OTF).
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Inputs.
ncpath_flat = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_36.nc"
ncpath_old  = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-old-001.nc"
ncpath_good = "/Users/zamperini/Documents/d3d_work/167196/167196-a2-tor240_44.nc"

# Modifier on the "old" data to bring it down some.
old_mod = 0.5

# If we want log yscales.
log = True

# Trim the tips of the probes where it shoots up erronously.
flat_start = 9
old_start  = 0
good_start = 9

def get_deps(ncpath):

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


a2_path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2_path)
a2_itf_x = a2["Distance from Tip D (cm)"].values
a2_otf_x = a2["Distance from Tip U (cm)"].values
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values

flat = get_deps(ncpath_flat)
good = get_deps(ncpath_good)
old  = get_deps(ncpath_old)

#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4))
fig = plt.figure(figsize=(11, 4))
gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)
(ax1, ax2, ax3) = gs.subplots(sharex='col', sharey='row')
ax11 = ax1.twinx()
ax22 = ax2.twinx()
ax33 = ax3.twinx()

ax11.scatter(a2_itf_x, a2_itf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:red")
ax11.scatter(a2_otf_x, a2_otf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple")
ax22.scatter(a2_itf_x, a2_itf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:red", label="ITF (RBS)")
ax22.scatter(a2_otf_x, a2_otf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple", label="OTF (RBS)")
ax33.scatter(a2_itf_x, a2_itf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:red", label="ITF (RBS)")
ax33.scatter(a2_otf_x, a2_otf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple", label="OTF (RBS)")

idx1 = old_start
idx2 = len(old["side1_x"]) - old_start - 1
ax1.plot(old["side1_x"][idx1:], old["side1_y"][idx1:]*old_mod, label="OTF", lw=2, color="tab:purple")
ax1.plot(old["side2_x"][:idx2], old["side2_y"][:idx2]*old_mod, label="ITF", lw=2, color="tab:red")

idx1 = good_start
idx2 = len(good["side1_x"]) - good_start - 1
ax2.plot(good["side1_x"][idx1:], good["side1_y"][idx1:], label="OTF (3DLIM)", lw=2, color="tab:purple")
ax2.plot(good["side2_x"][:idx2], good["side2_y"][:idx2], label="ITF (3DLIM)", lw=2, color="tab:red")

idx1 = flat_start
idx2 = len(flat["side1_x"]) - flat_start - 1
ax3.plot(flat["side1_x"][idx1:], flat["side1_y"][idx1:], label="OTF (3DLIM)", lw=2, color="tab:purple")
ax3.plot(flat["side2_x"][:idx2], flat["side2_y"][:idx2], label="ITF (3DLIM)", lw=2, color="tab:red")

#ax1.set_xlabel("Distance along probe (cm)", fontsize=12)
ax2.set_xlabel("Distance along probe (cm)", fontsize=12)
#ax3.set_xlabel("Distance along probe (cm)", fontsize=12)
ax33.set_ylabel(r"RBS W Areal Density (W/$\mathdefault{cm^2}$)", fontsize=12)
ax1.set_ylabel("3DLIM Deposition (arbitrary)", fontsize=12)

if log:
    pass
else:
    ax11.set_yticks([])
    ax22.set_yticks([])
    #ax2.set_yticks([])

legend_elements = [
  Line2D([0], [0], marker='*', color='k', label='ITF (RBS)', markerfacecolor='tab:red', markersize=15, lw=0),
  Line2D([0], [0], marker='*', color='k', label='OTF', markerfacecolor='tab:purple', markersize=15, lw=0),
  Line2D([0], [0], color="tab:red", label="ITF (3DLIM)", lw=2),
  Line2D([0], [0], color="tab:purple", label="OTF", lw=2)
  ]
ax33.legend(fontsize=12, handles=legend_elements)

ax1.text(0.15, 0.92, "a)", fontsize=12, transform=ax1.transAxes)
ax11.text(0.15, 0.92, "a)", fontsize=12, transform=ax11.transAxes)
ax2.text(0.15, 0.92, "b)", fontsize=12, transform=ax2.transAxes)
ax22.text(0.15, 0.92, "b)", fontsize=12, transform=ax22.transAxes)
ax3.text(0.15, 0.92, "c)", fontsize=12, transform=ax3.transAxes)
ax33.text(0.15, 0.92, "c)", fontsize=12, transform=ax33.transAxes)

if log:
    ax11.set_ylim([0.001, None])
    ax22.set_ylim([0.001, None])
    ax33.set_ylim([0.001, None])
    ax1.set_ylim([0.1, 40])
    ax2.set_ylim([0.1, 40])
    ax3.set_ylim([0.1, 40])
else:
    ax11.set_ylim([0, None])
    ax22.set_ylim([0, None])
    ax33.set_ylim([0, None])
    ax1.set_ylim([0, 40])
    ax2.set_ylim([0, 40])
    ax3.set_ylim([0, 40])

ax1.set_title("Rectangular (near-SOL accumulation)", fontsize=12, color="tab:pink")
ax2.set_title("DIVIMP (near-SOL accumulation)", fontsize=12, color="tab:green")
ax3.set_title("Flat (no near-SOL accumulation)", fontsize=12, color="tab:orange")

if log:
    ax1.set_yscale("log")
    ax11.set_yscale("log")
    ax2.set_yscale("log")
    ax22.set_yscale("log")
    ax3.set_yscale("log")
    ax33.set_yscale("log")

fig.tight_layout()
fig.show()
