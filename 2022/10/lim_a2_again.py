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

# Inputs. The first is from d3d-167196-blobby-002a.nc, while the second is from d3d-167196-mrc-shifted-nodrift-2.nc.
#ncpath_good = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240_46-blob.nc"
#ncpath_good = "/Users/zamperini/Documents/d3d_work/archive/167196/167196-a2-tor240_44.nc"
#ncpath_good = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240_44-blob-comp-001a.nc"
# ncpath_good = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-011.nc"
# ncpath_good = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-012.nc"
ncpath_good = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-013-018d.nc"

# THIS IS NOT A CORRECT WAY OF THINKING IT SEEMS.
exp_time = 4 * 25
absfac = 1e15
# absfac = 3e15
# absfac = 5.866e14
# absfac = 1.0
# absfac = 7.88999e17

# If we want log yscales.
log = False

# Trim the tips of the probes where it shoots up erronously.
good_start = 0

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


# Need to divide by exposure time to match 3DLIM units.
a2_path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2_path)
a2_itf_x = a2["Distance from Tip D (cm)"].values
a2_otf_x = a2["Distance from Tip U (cm)"].values
a2_itf_r = a2["R D (cm)"].values
a2_otf_r = a2["R U (cm)"].values
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values / exp_time
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values / exp_time

good = get_deps(ncpath_good)

fig, ax1 = plt.subplots(figsize=(5,4))
#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4))
#fig = plt.figure(figsize=(11, 4))
#gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)
#(ax1, ax2, ax3) = gs.subplots(sharex='col', sharey='row')
ax11 = ax1.twinx()


#ax11.scatter(a2_itf_x, a2_itf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:red")
#ax11.scatter(a2_otf_x, a2_otf_y, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple")
ax1.scatter(a2_itf_x, a2_itf_y*1e19, marker="*", edgecolors="k", s=200, zorder=5, color="tab:red")  # 1e15 W/cm2 to W/m2
ax1.scatter(a2_otf_x, a2_otf_y*1e19, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple")

idx1 = good_start
idx2 = len(good["side1_x"]) - good_start - 1
ax1.plot(good["side1_x"][idx1:], good["side1_y"][idx1:]*absfac, label="OTF (3DLIM)", lw=2, color="tab:purple")
ax1.plot(good["side2_x"][:idx2], good["side2_y"][:idx2]*absfac, label="ITF (3DLIM)", lw=2, color="tab:red")

ax1.set_ylabel("W Areal Density (W/m2/s)", fontsize=12)

if log:
    pass
else:
    ax11.set_yticks([])


legend_elements = [
  Line2D([0], [0], marker='*', color='k', label='ITF (RBS)', markerfacecolor='tab:red', markersize=15, lw=0),
  Line2D([0], [0], marker='*', color='k', label='OTF', markerfacecolor='tab:purple', markersize=15, lw=0),
  Line2D([0], [0], color="tab:red", label="ITF (3DLIM)", lw=2),
  Line2D([0], [0], color="tab:purple", label="OTF", lw=2)
  ]
ax11.legend(fontsize=12, handles=legend_elements)

ax1.text(0.15, 0.92, "a)", fontsize=12, transform=ax1.transAxes)
ax11.text(0.15, 0.92, "a)", fontsize=12, transform=ax11.transAxes)

# if log:
#     ax11.set_ylim([0.001, None])
#     ax1.set_ylim([0.1, 40])
# else:
#     ax11.set_ylim([0, None])
#     ax1.set_ylim([0, 40])


if log:
    ax1.set_yscale("log")
    ax11.set_yscale("log")

fig.tight_layout()
fig.show()
