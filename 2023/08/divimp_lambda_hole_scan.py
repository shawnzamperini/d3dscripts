# I did a scan in lambda_hole for DIVIMP to estimate how the lambda_hole value affects the core penetration.
import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pickle
import netCDF4
from scipy.interpolate import interp1d


root = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-020b-nocore"
cases = ["-l1", "-l2", "-l3", "", "-l4", "-l5", "-l6", "-l7", "-l8", "-l9", "-l10"]
ops = []
for case in cases:
    ops.append(oedge_plots.OedgePlots(root + "{}.nc".format(case)))

# The hole_lambda values.
hl = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.0150, 0.0175, 0.02, 0.03, 0.04, 0.05]

# Use the profile along the last core ring as our metric.
ircore = ops[0].irsep - 1
div_psins = []
div_nws = []
for i in range(0, len(ops)):
    op = ops[i]
    if i == 3:
        rings_only = True
    else:
        rings_only = False
    div_data = op.fake_probe(2.17, 2.36, -0.188, -0.188, "nz", charge="all", rings_only=rings_only)
    div_psin = np.array(div_data["psin"])
    nonan = ~np.isnan(div_psin)
    div_psin = div_psin[nonan]
    div_absfac = op.absfac
    div_nw = np.array(div_data["nz"])[nonan]
    mask = np.logical_and(div_psin >= 0.8, div_psin <= 1.2)
    div_psins.append(div_psin[mask])
    div_nws.append(div_nw[mask])

# Now load the SXR data. Copied from full_w_modeling, see comments there.
omfit_prof = netCDF4.Dataset("/Users/zamperini/Documents/d3d_work/divimp_files/167196/OMFITprofiles_167196_FIT.nc")
core_fit_psin = omfit_prof["psi_n"][:].data
core_fit_ne = omfit_prof["n_e"][:].mean(axis=0).data
core_fit_ne_err = np.sqrt(np.square(omfit_prof["n_e__uncertainty"][:].data).sum(axis=0))  # All zeros for some reason.
f_core_fit_ne = interp1d(core_fit_psin, core_fit_ne)
f_core_fit_ne_err = interp1d(core_fit_psin, core_fit_ne_err)
core_psin = np.linspace(0, 1.0, 300)
sxr_path = "/Users/zamperini/Documents/d3d_work/files/imp_analysis_167196.npz"
imps = np.load(sxr_path)
psin = np.square(imps["rho"])
wconc = imps["cw"]
core_mask = psin <= 0.2
wconc_avg = np.nanmean(wconc[:, core_mask])
wconc_std = np.nanstd(wconc[:, core_mask])
near_axis_psin = np.linspace(0.0, 0.2, 50)
nW_axis = np.mean(f_core_fit_ne(near_axis_psin) * wconc_avg)
nW_axis_min = np.mean(f_core_fit_ne(near_axis_psin) * (wconc_avg - wconc_std))
nW_axis_max = np.mean(f_core_fit_ne(near_axis_psin) * (wconc_avg + wconc_std))
tgyro_path = "/Users/zamperini/Documents/d3d_work/files/167196_tgyro_v1.pickle"
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

# One additional case comparing the results when parallel transport is ON.
op_par = oedge_plots.OedgePlots("/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-021.nc")
div_data = op_par.fake_probe(2.17, 2.36, -0.188, -0.188, "nz", charge="all", rings_only=rings_only)
div_psin = np.array(div_data["psin"])
nonan = ~np.isnan(div_psin)
div_psin = div_psin[nonan]
div_absfac = op_par.absfac
div_nw = np.array(div_data["nz"])[nonan]
mask = np.logical_and(div_psin >= 0.8, div_psin <= 1.2)
div_psin_par = div_psin[mask]
div_nw_par = div_nw[mask]

fig, ax = plt.subplots(figsize=(5, 4))
ax.axvline(1.0, color="k", linestyle="--")
ax.fill_between(tgyro_psin, tgyro_wdens_core_min, tgyro_wdens_core_max, color="tab:red", alpha=0.3)
ax.plot(tgyro_psin, tgyro_wdens_core, color="k", lw=3)
ax.plot(tgyro_psin, tgyro_wdens_core, color="tab:red", lw=2, label="SXR+TGYRO")
for i in range(0, len(ops)):

    # Make sure the OG case is on top of the others.
    if i == 3:
        ax.plot(div_psins[i], div_nws[i], color="k", lw=3, zorder=48)
        ax.plot(div_psins[i], div_nws[i], label=hl[i], color="tab:green", lw=2, zorder=49)
    else:
        ax.plot(div_psins[i], div_nws[i], label=hl[i], color="grey")
ax.plot(div_psin_par, div_nw_par, color="k", zorder=50, lw=3)
ax.plot(div_psin_par, div_nw_par, color="tab:pink", zorder=51, lw=2)
ax.set_yscale("log")
# ax.legend()
ax.set_xlim([0.75, 1.25])
ax.set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=14)
ax.set_ylabel(r"W Density ($\mathdefault{m^{-3}}$)", fontsize=14)
ax.grid()
ax.arrow(0.92, 2e12, 0.0, 7e14, width=0.005, color="k", zorder=100, head_length=5e14)
ax.text(0.82, 1e12, r"$\mathdefault{\lambda_{hole}}$ = 0.25 cm", backgroundcolor="w", fontsize=12)
ax.text(0.90, 2e15, "5 cm", backgroundcolor="w", fontsize=12)
ax.text(1.04, 2e15, "|| transport in blob", fontsize=12, backgroundcolor="w")
ax.text(1.07, 6e14, "ON", color="tab:pink", fontsize=12, backgroundcolor="w")
ax.text(1.14, 6e14, "OFF", color="tab:green", fontsize=12, backgroundcolor="w")
fig.tight_layout()
fig.show()
