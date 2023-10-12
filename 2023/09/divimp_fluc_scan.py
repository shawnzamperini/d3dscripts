import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
import pickle
import netCDF4
from scipy.interpolate import interp1d, griddata
import LimPlots


# Load all the cases.
root = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/"
epolscan = [100, 500, 1500, 2000]
# probscan = [0.1, 0.25, 0.75, 0.9, 0.4, 0.45, 0.55, 0.6]
probscan = [0.4, 0.45, 0.55, 0.6]
problabels = [5, 6, 7, 8]
ops_epol = []
ops_prob = []
for i in range(1, len(epolscan) + 1):
    print("{}/{}".format(i, len(epolscan)))
    case1 = "d3d-167196-fluc-002-epolscan{}.nc".format(i)
    ops_epol.append(oedge_plots.OedgePlots(root + case1))
for i in range(0, len(probscan)):
    print("{}/{}".format(i, len(probscan)))
    case2 = "d3d-167196-fluc-002-probscan{}.nc".format(problabels[i])
    ops_prob.append(oedge_plots.OedgePlots(root + case2))
op_base = oedge_plots.OedgePlots(root + "d3d-167196-fluc-002.nc")

# Load the SXR + TGYRO results. Just big hunk of code here, removed the comments from full_w_modeling.py.
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
tgyro_path = "/Users/zamperini/Documents/d3d_work/files/167196_tgyro_v1.pickle"
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

# Likewise for 3DLIM data.
lim_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-013-018d-noprobe.nc"
lp = LimPlots.LimPlots(lim_path)
lpdata = lp.plot_par_rad("nz", 21, charge="all")
mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 15
absfac = 1e15
print("Creating psin interpolation...")
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
rorigin = 2.282
rad_locs = rorigin - lpdata["Y"][0]
lim_nz = lpdata["Z"][mididx].data * absfac
lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
lim_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(lim_rzs))
wall_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), (rorigin + 0.08, -0.188))
mask_lim = lim_psins >= 1.18
lim_psins = lim_psins[mask_lim]
lim_nz = lim_nz[mask_lim]


# Function to get radial W profiles.
def get_rad_w_prof(op):
    div_data = op.fake_probe(2.17, 2.36, -0.188, -0.188, "nz", charge="all", rings_only=False)
    div_psin = np.array(div_data["psin"])
    nonan = ~np.isnan(div_psin)
    div_psin = div_psin[nonan]
    div_absfac = op.absfac
    div_nw = np.array(div_data["nz"])[nonan]
    mask = np.logical_and(div_psin >= 0.8, div_psin <= 1.2)
    div_psin = div_psin[mask]
    div_nw = div_nw[mask]
    return div_psin, div_nw


# Now the radial profiles from the DIVIMP runs.
psins_epol = []
nws_epol = []
psins_prob = []
nws_prob = []
for i in range(0, len(epolscan)):
    psin, nw = get_rad_w_prof(ops_epol[i])
    psins_epol.append(psin)
    nws_epol.append(nw)
for i in range(0, len(probscan)):
    psin, nw = get_rad_w_prof(ops_prob[i])
    psins_prob.append(psin)
    nws_prob.append(nw)
psin_base, nw_base = get_rad_w_prof(op_base)


def label_at_psin(ax, x, y, label, psin):
    idx = np.argmin(np.abs(x - psin))
    ax.text(x[idx], y[idx], label, fontsize=7, bbox={"facecolor":"w", "edgecolor":"k"}, horizontalalignment="center")


# Plots assembling everything.
psin_labels = [1.07, 1.13, 1.15, 1.15, 1.15]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
ax1.set_xlim([0.8, 1.3])
for ax in [ax1, ax2]:
    ax.axvline(1.0, linestyle="--", color="k")
    ax.axvline(1.4, color="k")
    ax.fill_between(tgyro_psin, tgyro_wdens_core_min, tgyro_wdens_core_max, color="tab:red", alpha=0.3)
    ax.plot(tgyro_psin, tgyro_wdens_core, color="k", lw=3)
    ax.plot(tgyro_psin, tgyro_wdens_core, color="tab:red", lw=2, label="SXR+TGYRO")
    ax.plot(lim_psins, lim_nz, color="k", lw=3)
    ax.plot(lim_psins, lim_nz, color="tab:green", lw=2, label="3DLIM")
    ax.set_yscale("log")
    ax.set_ylim([1e9, 1e16])
    ax.grid(alpha=0.3)
    ax.set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=14)
for i in range(0, len(epolscan)):
    ax1.plot(psins_epol[i], nws_epol[i], lw=2, linestyle="-", color="tab:pink")
for i in range(0, len(probscan)):
    ax2.plot(psins_prob[i], nws_prob[i], lw=2, linestyle="-", color="tab:pink")
    label_at_psin(ax2, psins_prob[i], nws_prob[i], probscan[i], psin_labels[i])
ax1.plot(psin_base, nw_base, color="k", lw=3)
ax1.plot(psin_base, nw_base, color="tab:pink", lw=2)
ax2.plot(psin_base, nw_base, color="k", lw=3)
ax2.plot(psin_base, nw_base, color="tab:pink", lw=2)
label_at_psin(ax2, psin_base, nw_base, 0.5, 1.15)
# ax1.arrow(1.05, 1e11, (1.18-1.05), (5e14 - 1e11), width=0.005, head_length=1e14, zorder=100, color="k")
ax1.text(1.05, 1e15, r"$\mathdefault{E_{\theta}}}$ = 100-2,000 V/m", fontsize=10)
ax1.text(0.81, 4e15, "a)", fontsize=10)
ax2.text(0.81, 4e15, "b)", fontsize=10)
ax1.set_ylabel(r"W Density ($\mathdefault{m^{-3}}$)", fontsize=14)
fig.tight_layout()
fig.show()
