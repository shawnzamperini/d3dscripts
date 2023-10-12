import oedge_plots
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
import LimPlots
from scipy.interpolate import griddata
import pickle

plt.rcParams["font.family"] = "Century Gothic"
plt.rc('axes', unicode_minus=False)

# Selecting a decent 3DLIM run and the needed ABSFAC.
limpath = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-011.nc"
exp_time = 4 * 25
absfac = 1e17 / exp_time  # Showing that we needed to divide by the exposure time so RBS and 3DLIM units agree.

# Selecting the DIVIMP runs.
diffpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-diff-001-predep.nc"
# blobpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010f-predep.nc"
blobpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010f-old-bkg.nc"
diff = oedge_plots.OedgePlots(diffpath)
blob = oedge_plots.OedgePlots(blobpath)

# If we want log yscales.
log = False


def get_deps(ncpath):
    # Load netcdf file.
    nc = netCDF4.Dataset(ncpath)

    # Location of each P bin, and its width. Note syntax in pulling out data.
    ps = nc.variables['PS'][:].data
    pwids = nc.variables['PWIDS'][:].data

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

    return {"side1_x": side1_x, "side1_y": side1_y, "side2_x": side2_x, "side2_y": side2_y}


# Load RBS data.
a2_path = "/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx"
a2 = pd.read_excel(a2_path)
a2_itf_x = a2["Distance from Tip D (cm)"].values
a2_otf_x = a2["Distance from Tip U (cm)"].values
a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values
a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values

lim = get_deps(limpath)

# Load the gfile and map the R, Z points to psin.
print("Creating psin interpolation...")
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]

def get_lim_data(lim_path, lim_absfac):
    lp = LimPlots.LimPlots(lim_path)
    lpdata = lp.plot_par_rad("nz", 21, charge="all")

    # Let's not choose right in the middle (Y = 0) since flows are zero there, and it actually leads to distracting
    # valleys in the radial profile due to no impurities there.
    mididx = np.argmin(np.abs(lpdata["X"][:, 0])) + 5

    rorigin = 2.282
    rad_locs = rorigin - lpdata["Y"][0]

    # This value is chosen in lim_a2_again.py. It is that needed to match the experimental deposition.
    lim_nz = lpdata["Z"][mididx].data * lim_absfac
    lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
    lim_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(lim_rzs))
    lim_rhos = np.sqrt(lim_psins)
    wall_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), (rorigin + 0.08, -0.188))
    wall_rho = np.sqrt(wall_psin)
    return {"psin":lim_psins, "nz":lim_nz}

lim_blob = get_lim_data("/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-011-noprobe.nc", absfac)

# Similar to how I did in tomas_sxr_results.py, I will compare the radial profiles from DIVIMP vs the one from 3DLIM
# to see how they line up in magnitude. I will compare the normal vs. the blobby DIVIMP ones above to a 3DLIM run where
# the necessary absfac has been determined based off what gives the best match to the profiles in magnitude.
diff_probe = diff.fake_probe(2.21, 2.36, -0.188, -0.188, "nz", charge="all")
blob_probe = blob.fake_probe(2.21, 2.36, -0.188, -0.188, "nz", charge="all")

div_mask = np.array(blob_probe["psin"]) < 1.14
# div_mask = np.array(blob_probe["psin"]) < 99
div_psin = np.array(blob_probe["psin"])[div_mask]
div_blob_nz = np.array(blob_probe["nz"])[div_mask]
div_diff_nz = np.array(diff_probe["nz"])[div_mask]


# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
fig, ax1 = plt.subplots(figsize=(5, 4))
ax11 = ax1.twinx()

ax1.scatter(a2_itf_x[:-1], a2_itf_y[:-1] * 1e19, marker="*", edgecolors="k", s=200, zorder=5,
            color="tab:red")  # 1e15 W/cm2 to W/m2
ax1.scatter(a2_otf_x[:-3], a2_otf_y[:-3] * 1e19, marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple")

# Trim the tips of the probes where it shoots up erronously.
good_start = 6
idx1 = good_start
idx2 = len(lim["side1_x"]) - good_start - 1
ax1.plot(lim["side1_x"][idx1:], lim["side1_y"][idx1:] * absfac * exp_time, label="OTF (3DLIM)", lw=2, color="tab:purple")
ax1.plot(lim["side2_x"][:idx2], lim["side2_y"][:idx2] * absfac * exp_time, label="ITF (3DLIM)", lw=2, color="tab:red")
ax1.set_xlabel("Distance along probe (cm)", fontsize=12)
ax1.set_ylabel(r"W Areal Density $\mathdefault{(W/m^2)}$", fontsize=12)
ax1.set_ylim([0, 5.5e18])

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

# ax11.text(0.02, 0.94, "a)", fontsize=12, transform=ax11.transAxes)

if log:
    ax1.set_yscale("log")
    ax11.set_yscale("log")

fig.tight_layout()
fig.show()

fig, ax2 = plt.subplots(figsize=(5, 4))

ax2.axvline(1.0, color="k")
ax2.plot(div_psin, div_blob_nz, label="Blobby", color="tab:red", lw=3)
ax2.plot(div_psin, div_diff_nz, label="Diffusive", color="tab:purple", lw=3)
ax2.plot(lim_blob["psin"], lim_blob["nz"], label="3DLIM", color="tab:cyan", lw=3)
ax2.set_yscale("log")
ax2.legend(fontsize=12)
ax2.grid(alpha=0.3)
ax2.set_xlabel(r"$\mathdefault{\psi_N}$", fontsize=12)
ax2.set_ylabel(r"W Density $\mathdefault{(m^{-3})}$", fontsize=12)
# ax2.text(0.02, 0.94, "b)", fontsize=12, transform=ax2.transAxes)

fig.tight_layout()
fig.show()
