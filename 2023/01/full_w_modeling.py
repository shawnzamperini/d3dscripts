# A new script to combine experimental and modeling results for a full-SOL picture of the W density.
import pickle
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.signal import savgol_filter, medfilt
import pandas as pd
import oedge_plots
import LimPlots
from matplotlib import lines

# Flags to turn on different plots.
plot1 = False
plot2 = True
plot3 = True
shot = 167196 # Either 167196 or 190423

# First load fit profiles of the electron density. We will use this in our assumption of "The W concentration stays
# constant throughout the core". The uncertainties of the average values requires some error propagation.
if shot == 167196:
    omfit_prof = netCDF4.Dataset("/Users/zamperini/Documents/d3d_work/divimp_files/167196/OMFITprofiles_167196_FIT.nc")
elif shot == 190423:
    omfit_prof = netCDF4.Dataset("/Users/zamperini/Documents/d3d_work/divimp_files/190423/OMFITprofiles_190423_FIT.nc")
core_fit_psin = omfit_prof["psi_n"][:].data
core_fit_ne = omfit_prof["n_e"][:].mean(axis=0).data
core_fit_ne_err = np.sqrt(np.square(omfit_prof["n_e__uncertainty"][:].data).sum(axis=0))  # All zeros for some reason.

# Interpolations of the density for the below SXR analysis.
f_core_fit_ne = interp1d(core_fit_psin, core_fit_ne)
f_core_fit_ne_err = interp1d(core_fit_psin, core_fit_ne_err)
core_psin = np.linspace(0, 1.0, 300)

# Now load the SXR analysis provided by Tomas.
print("Loading SXR data...")
if shot == 167196:
    sxr_path = "/Users/zamperini/Documents/d3d_work/files/imp_analysis_167196.npz"
    imps = np.load(sxr_path)
    psin = np.square(imps["rho"])
    wconc = imps["cw"]
    core_mask = psin <= 0.2
    wconc_avg = np.nanmean(wconc[:, core_mask])
    wconc_std = np.nanstd(wconc[:, core_mask])

    tgyro_path = "/Users/zamperini/Documents/d3d_work/files/167196_tgyro_v1.pickle"


elif shot == 190423:
    with open("/Users/zamperini/Documents/d3d_work/files/190422_w_imprad.pickle", "rb") as f:
        sxr = pickle.load(f)
    sxr_t = sxr.t.data
    sxr_psin = np.square(sxr.rho.data)
    sxr_wdens = sxr["W dens."].data

    # Restrict the SXR data between a specified time range that can considered constant. Inspecting the data shows that
    # a range of 2000-4500 is a good range.
    sxr_mask = np.logical_and(sxr_t >= 2000, sxr_t <= 4500)
    sxr_t = sxr_t[sxr_mask]
    tmp = sxr_wdens[sxr_mask].mean(axis=0)  # I think this will propagate the uncertainty values.
    sxr_wdens = np.array([v.n for v in tmp])
    sxr_wdens_err = np.array([v.s for v in tmp])

    # Create interpolation functions of ne(psin) and wdens(psin) so that we can plot the ratio of them on a common psin.
    f_sxr_wdens = interp1d(sxr_psin, sxr_wdens)
    f_sxr_wdens_err = interp1d(sxr_psin, sxr_wdens_err)
    wconc = f_sxr_wdens(core_psin) / f_core_fit_ne(core_psin)

    # Just use the values in the deep core as representative. Define here as psin = [0.0, 0.2].
    core_mask = core_psin <= 0.2
    wconc_avg = (f_sxr_wdens(core_psin[core_mask]) / f_core_fit_ne(core_psin[core_mask])).mean()
    wconc_std = (f_sxr_wdens(core_psin[core_mask]) / f_core_fit_ne(core_psin[core_mask])).std()

    tgyro_path = "/Users/zamperini/Documents/d3d_work/files/190423_tgyro_v1.pickle"

near_axis_psin = np.linspace(0.0, 0.2, 50)
nW_axis = np.mean(f_core_fit_ne(near_axis_psin) * wconc_avg)
nW_axis_min = np.mean(f_core_fit_ne(near_axis_psin) * (wconc_avg - wconc_std))
nW_axis_max = np.mean(f_core_fit_ne(near_axis_psin) * (wconc_avg + wconc_std))

with open(tgyro_path, "rb") as f:
    tgyro = pickle.load(f)
n_blend_i1 = tgyro["n_blend_i1"]
n_ratio = tgyro["n_ratio"]
tgyro_rho = tgyro["rho"]
tgyro_psin = np.square(tgyro_rho)
nW0 = n_blend_i1 / n_ratio
tgyro_mask = tgyro_psin <= 0.2
tgyro_wdens_core = nW0 / nW0[tgyro_mask].mean() * nW_axis
tgyro_wdens_core_min = nW0 / nW0[tgyro_mask].mean() * nW_axis_min
tgyro_wdens_core_max = nW0 / nW0[tgyro_mask].mean() * nW_axis_max

# Then get the estimated radial profile of W in the core.
print("Average W concentration: {:.2e} +/- {:.2e}".format(wconc_avg, wconc_std))
wdens_core = f_core_fit_ne(core_psin) * wconc_avg
wdens_core_min = f_core_fit_ne(core_psin) * (wconc_avg - wconc_std)
wdens_core_max = f_core_fit_ne(core_psin) * (wconc_avg + wconc_std)

if plot1 and shot == 190423:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

    ax1.axvline(1.0, color="k")
    ax1.fill_between(core_psin, f_sxr_wdens(core_psin) - f_sxr_wdens_err(core_psin),
                     f_sxr_wdens(core_psin) + f_sxr_wdens_err(core_psin), color="tab:pink", alpha=0.3)
    ax1.plot(core_psin, f_sxr_wdens(core_psin), label="SXR", color="tab:pink", lw=3)
    ax11 = ax1.twinx()
    ax11.fill_between(core_psin, f_core_fit_ne(core_psin) - f_core_fit_ne_err(core_psin),
                      f_core_fit_ne(core_psin) + f_core_fit_ne_err(core_psin), color="tab:cyan", alpha=0.3)
    ax11.plot(core_psin, f_core_fit_ne(core_psin), label="ne", color="tab:cyan", lw=3)
    ax1.set_ylim([0, 1e17])
    ax1.set_xlabel("Psin", fontsize=14)
    ax1.set_ylabel("W Density (m-3)", fontsize=14, color="tab:pink")
    ax11.set_ylabel("ne (m-3)", fontsize=14, color="tab:cyan")

    ax2.axvline(1.0, color="k")
    ax2.plot(core_psin, wconc, color="tab:red", lw=3)
    ax2.axhline(wconc_avg, color="tab:red", linestyle="--")
    ax2.set_ylabel("W Concentration", fontsize=14)
    ax2.set_xlabel("Psin", fontsize=14)
    ax2.set_ylim([1e-5, 1e-3])
    ax2.grid(which="both", alpha=0.3)
    ax2.set_yscale("log")

    fig.tight_layout()
    fig.show()


# Now load a 3DLIM run and get the centerline data.
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


# Extract the arrays, assign accordingly.
if shot == 167196:
    absfac = 1e15
    exposure_time = 4 * 25  # 4 second flattop, 25 repeat shots.
    # deps_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-011.nc"  # Case with blobby-010f source.
    deps_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-013-018d.nc" # Case with blobby-018d source.
elif shot == 190423:
    absfac = 2.2e15  # In units of W / s
    exposure_time = 4 * 2
    # deps_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-004.nc"
    deps_path = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-mcp01w-002.nc"
deps = get_deps(deps_path)
lim_itfx = deps["side2_x"]
lim_itfy = deps["side2_y"] * absfac * exposure_time
lim_otfx = deps["side1_x"]
lim_otfy = deps["side1_y"] * absfac * exposure_time
window = 21
lim_itfy = savgol_filter(medfilt(lim_itfy, window), window, 2)
lim_otfy = savgol_filter(medfilt(lim_otfy, window), window, 2)

if shot == 167196:
    itf_mask = lim_itfx > 0.5
    otf_mask = lim_otfx > 0.5
elif shot == 190423:
    itf_mask = lim_itfx > 0.0
    otf_mask = lim_otfx > 0.0
lim_itfx = lim_itfx[itf_mask]
lim_itfy = lim_itfy[itf_mask]
lim_otfx = lim_otfx[otf_mask]
lim_otfy = lim_otfy[otf_mask]

# Load the RBS data. The y data NEEDS TO BE DIVIDED BY EXPOSURE TIME FOR 3DLIM COMPARISON.
if shot == 167196:
    a2 = pd.read_excel("/Users/zamperini/My Drive/School/Tennessee/Research/Collector Probe Excel Sheets/A2.xlsx")
    a2_itf_x = a2["Distance from Tip D (cm)"].values
    a2_otf_x = a2["Distance from Tip U (cm)"].values
    a2_itf_y = a2["W Areal Density D (1e15 W/cm2)"].values * 1e19  # In units of W/m2/s
    a2_otf_y = a2["W Areal Density U (1e15 W/cm2)"].values * 1e19

# The LAMS data for MCP01W.
elif shot == 190423:
    itf_path = "/Users/zamperini/My Drive/Research/Data/cp_data/MCPL01W.csv"
    otf_path = "/Users/zamperini/My Drive/Research/Data/cp_data/MCPR01W.csv"
    itf = pd.read_csv(itf_path)
    otf = pd.read_csv(otf_path)
    lam_itfx = itf["Distance along probe (cm)"].values
    lam_itfy = itf["W areal density (1e15cm-2)"].values * 1e15 * 10000
    lam_otfx = otf["Distance along probe (cm)"].values
    lam_otfy = otf["W areal density (1e15cm-2)"].values * 1e15 * 10000
    window = 21
    lam_itfy = savgol_filter(medfilt(lam_itfy, window), window, 2)
    lam_otfy = savgol_filter(medfilt(lam_otfy, window), window, 2)

if plot2:
    fig, ax = plt.subplots(figsize=(5, 4))
    if shot == 167196:
        ax.scatter(a2_itf_x[:-1], a2_itf_y[:-1], marker="*", edgecolors="k", s=200, zorder=5,
                    color="tab:red", label="ITF (RBS)")
        ax.scatter(a2_otf_x[:-3], a2_otf_y[:-3], marker="*", edgecolors="k", s=200, zorder=5, color="tab:purple",
                   label="OTF")
        ax.set_ylim([0, 5.5e16])
    elif shot == 190423:
        ax.plot(lam_itfx, lam_itfy, color="tab:red", lw=3, label="LAMS")
        ax.plot(lam_otfx, lam_otfy, color="tab:purple", lw=3)
        ax.set_ylim([0, 2e17])
    ax.plot(lim_itfx, lim_itfy, color="k", lw=5)
    ax.plot(lim_otfx, lim_otfy, color="k", lw=5)
    ax.plot(lim_itfx, lim_itfy, label="ITF (3DLIM)", color="tab:red", lw=3)
    ax.plot(lim_otfx, lim_otfy, color="tab:purple", lw=3, label="OTF")
    # ax.legend(font-size=14)
    ax.set_xlabel("Distance along probe (cm)", fontsize=14)
    ax.set_ylabel(r"W Areal Density $\mathdefault{(W/m^2)}$", fontsize=14)
    # ax.set_title("ABSFAC = {:.2e}".format(absfac))
    if shot == 167196:
        ax.set_xlim([0, 7])
        ax.set_ylim([0, 5.5e18])
    elif shot == 190423:
        ax.set_xlim([0, 6])
    ax.legend(fontsize=12)

    fig.tight_layout()
    fig.show()

# Load a DIVIMP case so that we can overlay that on our master plot as well.
if shot == 167196:
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-021.nc"
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-020b-nocore.nc"
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-018d.nc"
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-010f.nc"
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-013f-0.4exb.nc"
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-fluc-002.nc"  # tcorr = 1.5e-6
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-fluc-002-tcorrscan1.nc"  # tcorr = 3e-6, this is just simpler than 1.5, cleaner, otherwise the same results

    op = oedge_plots.OedgePlots(ncpath)
    div_data = op.fake_probe(2.17, 2.36, -0.188, -0.188, "nz", charge="all", rings_only=False)

elif shot == 190423:
    # tungsten-002: ring-entry mach flow, 0.3 m2/s + blobby @ 1.5 kHz
    # tungsten-003: ring-entry mach flow, 0.3 m2/s + blobby @ 3.0 kHz
    # tungsten-004: ring-entry mach flow, 0.03 m2/s + blobby @ 3.0 kHz
    # tunsgten-005: ring-entry mach flow, 0.00 m2/s + blobby from 167196-blobby-004 @ 2 kHz (Gaussian blob vr at 800 m/s, 400 m/s width)
    # tungsten-006: same as 005, but drifts ON.
    # tungsten-008: blobby model revisted after upgrades
    # tungsten-009: messing with the sputtering yield multiplier and corr time for less steep profile compared to 008.
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-tungsten-009.nc"
    op = oedge_plots.OedgePlots(ncpath)
    div_data = op.fake_probe(2.20, 2.30, -0.188, -0.188, "nz", charge="all", rings_only=True)

div_psin = np.array(div_data["psin"])
nonan = ~np.isnan(div_psin)
div_psin = div_psin[nonan]
div_absfac = op.absfac
div_nw = np.array(div_data["nz"])[nonan]
if shot == 167196:
    mask = np.logical_and(div_psin >= 0.8, div_psin <= 1.2)
elif shot == 190423:
    mask = np.logical_and(div_psin >= 0.8, div_psin <= 1.09)
div_psin = div_psin[mask]
div_nw = div_nw[mask]

# An additional comparison run without the blobby contribution.
if shot == 167196:
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-diff-001-predep-test.nc"
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-diff-005-predep.nc"
    # ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d-167196-modE-expW-dft60-W_test.nc"
    op = oedge_plots.OedgePlots(ncpath)
    div_data_diff = op.fake_probe(2.17, 2.36, -0.188, -0.188, "nz", charge="all", rings_only=False)
elif shot == 190423:
    # tungsten-001: ring-entry mach flow, 0.3 m2/s
    # tungsten-008: ring-entry mach flow (incl. core), 0.3 m2/2
    ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/d3d-190423-tungsten-008-diff.nc"
    op = oedge_plots.OedgePlots(ncpath)
    div_data_diff = op.fake_probe(2.20, 2.30, -0.188, -0.188, "nz", charge="all")

div_psin_diff = np.array(div_data_diff["psin"])
nonan_diff = ~np.isnan(div_psin_diff)
div_psin_diff = div_psin_diff[nonan_diff]
absfac_diff = op.absfac
div_nw_diff = np.array(div_data_diff["nz"])[nonan_diff]
if shot == 167196:
    mask_diff = np.logical_and(div_psin_diff >= 0.8, div_psin_diff <= 1.2)
elif shot == 190423:
    mask_diff = np.logical_and(div_psin_diff >= 0.8, div_psin_diff <= 1.2)
div_psin_diff = div_psin_diff[mask_diff]
div_nw_diff = div_nw_diff[mask_diff]

# Load a 3DLIM run, pull some radial profiles from different locations to
# capture 3D effects.
if shot == 167196:
    # lim_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-011-noprobe.nc"
    lim_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-013-018d-noprobe.nc"
elif shot == 190423:
    lim_path = "/Users/zamperini/Documents/d3d_work/lim_runs/190423/190423-noprobe-002.nc"
lp = LimPlots.LimPlots(lim_path)
lpdata = lp.plot_par_rad("nz", 21, charge="all")

# Let's not choose right in the middle (Y = 0) since flows are zero there, and it actually leads to distracting
# valleys in the radial profile due to no impurities there.
mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 15

# Need the following to map 3DLIM onto psin.
print("Creating psin interpolation...")
if shot == 167196:
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
elif shot == 190423:
    gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190423/190423_3000.pickle"

with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]

if shot == 167196:
    rorigin = 2.282
elif shot == 190423:
    rorigin = 2.305
rad_locs = rorigin - lpdata["Y"][0]
lim_nz = lpdata["Z"][mididx].data * absfac
lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
lim_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(lim_rzs))
wall_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), (rorigin + 0.08, -0.188))

# Ignore tip where we're digging into the source, not what we want.
if shot == 167196:
    mask_lim = lim_psins >= 1.18
elif shot == 190423:
    pass # Need to figure out still.
    mask_lim = lim_psins >= 0.0
lim_psins = lim_psins[mask_lim]
lim_nz = lim_nz[mask_lim]

# Load the MAFOT data of the connection length at the CP location.
if shot == 167196:
    mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/lam_mimes_plunge.dat"
elif shot == 190423:
    mafot_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190423/lam_cp_conns.dat"
def load_mafot(path):
    print("Loading MAFOT data...")
    columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
               "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
    mafot = pd.read_csv(path, skiprows=52, names=columns, delimiter="\t")
    #mafot_otf = pd.read_csv(otf_path, skiprows=52, names=columns, delimiter="\t")
    conns_r = mafot["R (m)"]
    conns_psin = mafot["psi"]
    conns = mafot["Lconn (km)"].values * 1000  # km to m
    #conns_l_otf = mafot_otf["Lconn (km)"].values * 1000
    return {"r":conns_r, "conns":conns, "psin":conns_psin}
mafot = load_mafot(mafot_path)

smooth = False
if plot3:

    fig, ax1 = plt.subplots(figsize=(7, 4))

    ax1.axvline(1.0, linestyle="--", color="k")
    ax1.axvline(1.4, color="k")
    #if shot == 167196:
    ax1.fill_between(tgyro_psin, tgyro_wdens_core_min, tgyro_wdens_core_max, color="tab:red", alpha=0.3)
    ax1.plot(tgyro_psin, tgyro_wdens_core, color="k", lw=3)
    ax1.plot(tgyro_psin, tgyro_wdens_core, color="tab:red", lw=2, label="SXR+TGYRO")
    # ax1.plot(core_psin, wdens_core, color="tab:red", lw=3, linestyle="--", label=r"$\mathdefault{n_W=c_Wn_e}$")
    #elif shot == 190423:
        #ax1.fill_between(core_psin, wdens_core_min, wdens_core_max, color="tab:red", alpha=0.3)
    #    ax1.plot(core_psin, wdens_core, color="tab:red", lw=3, linestyle="--", label=r"$\mathdefault{n_W=c_Wn_e}$")
    if smooth:
        ax1.plot(div_psin, savgol_filter(div_nw, 21, 2), color="k", lw=3)
        ax1.plot(div_psin, savgol_filter(div_nw, 21, 2), color="tab:pink", lw=2, label="DIVIMP (fluctuation)")
        ax1.plot(div_psin_diff, savgol_filter(div_nw_diff, 21, 2), color="k", lw=3)
        ax1.plot(div_psin_diff, savgol_filter(div_nw_diff, 21, 2), color="tab:pink", lw=2, linestyle="--")
    else:
        ax1.plot(div_psin, div_nw, color="k", lw=3)
        ax1.plot(div_psin, div_nw, color="tab:pink", lw=2, label="DIVIMP (fluctuation)")
        # ax1.fill_between(div_psin, div_nw*0.8, div_nw*1.2, color="tab:pink", label="DIVIMP")
        ax1.plot(div_psin_diff, div_nw_diff, color="k", lw=3)
        ax1.plot(div_psin_diff, div_nw_diff, color="tab:pink", lw=2, linestyle="--", label="DIVIMP (diffusive)")
        # ax1.fill_between(div_psin_diff, div_nw_diff*0.8, div_nw_diff*1.2, color="tab:cyan")
        pass
    ax1.plot(lim_psins, lim_nz, color="k", lw=3)
    ax1.plot(lim_psins, lim_nz, color="tab:green", lw=2, label="3DLIM")
    ax1.set_yscale("log")
    # ax1.legend(fontsize=14)
    # ax1.set_xlim([0.95, 1.45])
    ax1.set_xlim([0, 1.45])
    if shot == 167196:
        ax1.set_ylim([1e9, 1e16])
    elif shot == 190423:
        ax1.set_ylim(1e10, 5e16)
    ax1.set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=14)
    ax1.set_ylabel(r"W Density ($\mathdefault{m^{-3}}$)", fontsize=14)
    ax1.grid()

    # The legend.
    dotted_line1 = lines.Line2D([], [], linewidth=2, linestyle="--", color="tab:pink")
    dotted_line2 = lines.Line2D([], [], linewidth=3, linestyle="-", color="k")
    line1a = lines.Line2D([], [], linewidth=3, color="k")
    line1b = lines.Line2D([], [], linewidth=2, color="tab:pink")
    line2a = lines.Line2D([], [], linewidth=3, color="k")
    line2b = lines.Line2D([], [], linewidth=2, color="tab:red")
    line3a = lines.Line2D([], [], linewidth=3, color="k")
    line3b = lines.Line2D([], [], linewidth=2, color="tab:green")
    ax1.legend([(line2a, line2b), (line1a, line1b), (dotted_line2, dotted_line1), (line3a, line3b)],
               ["SXR+TGYRO", "DIVIMP (fluctuation)", "DIVIMP (diffusive)", "3DLIM"], fontsize=14)

    #ax11 = ax1.twinx()
    #ax11.plot(mafot["psin"], mafot["conns"], color="tab:purple", lw=3)
    #ax11.set_yscale("log")

    fig.tight_layout()
    fig.show()
