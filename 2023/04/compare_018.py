import oedge_plots
import matplotlib.pyplot as plt
import pickle
import numpy as np
import netCDF4
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


# First load fit profiles of the electron density. We will use this in our assumption of "The W concentration stays
# constant throughout the core". The uncertainties of the average values requires some error propagation.
omfit_prof = netCDF4.Dataset("/Users/zamperini/Documents/d3d_work/divimp_files/167196/OMFITprofiles_167196_FIT.nc")
core_fit_psin = omfit_prof["psi_n"][:].data
core_fit_ne = omfit_prof["n_e"][:].mean(axis=0).data
core_fit_ne_err = np.sqrt(np.square(omfit_prof["n_e__uncertainty"][:].data).sum(axis=0))  # All zeros for some reason.

# Interpolations of the density for the below SXR analysis.
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

div_psins = []; div_nws = []
case = "020"
if case == "018":
    labels = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
elif case in ["019"]:
    labels = ["a", "b", "c"]
elif case == "020":
    labels = ["a", "b", "c", "d"]
for label in labels:
    print(label)
    if case == "020":
        path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-{}{}-nocore.nc".format(case, label)
    else:
        path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-blobby-{}{}.nc".format(case, label)
    op = oedge_plots.OedgePlots(path)
    div_data = op.fake_probe(2.17, 2.36, -0.188, -0.188, "nz", charge="all", rings_only=True)
    div_psins.append(div_data["psin"])
    div_nws.append(div_data["nz"])


fig, ax = plt.subplots()
ax.axvline(1.0, color="k", linestyle="--")

ax.plot(tgyro_psin, tgyro_wdens_core, color="k", lw=3)
ax.plot(tgyro_psin, tgyro_wdens_core, color="tab:red", lw=2, label="SXR+TGYRO")

smooth = False
window = 15
for i in range(0, len(div_psins)):
    if smooth:
        ax.plot(div_psins[i][:-window], savgol_filter(div_nws[i], window, 2)[:-window], label=labels[i])
    else:
        ax.plot(div_psins[i], div_nws[i], label=labels[i])

ax.grid(alpha=0.3, which="both")
ax.legend()
ax.set_xlim([0.8, 1.4])
ax.set_yscale("log")
ax.set_xlabel("Psin")
ax.set_ylabel("nW (m-3)")
fig.tight_layout()
fig.show()