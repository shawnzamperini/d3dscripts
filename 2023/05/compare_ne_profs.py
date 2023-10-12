# Compare SOL density profiles from 167196 to 190422/3 (using 167195 and 190411 as a proxy).
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d
from scipy.optimize import curve_fit
import numpy as np
import pickle
import netCDF4


def load_data(shot):
    # Shot-specific constants.
    if shot == 167195:
        rcp_shift = -0.018
        print("RCP shifted by {:.3} m".format(rcp_shift))
        rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP167195_2.tab"
        gfile_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/167195_4460.pickle"
        omfit_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/OMFITprofiles_167196_FIT.nc"
    elif shot == 190423:
        rcp_shift = 0.0
        rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190411_2.tab"  # Assuming this data is good enough for 190423.
        gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190423/190423_3000.pickle"
        omfit_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/OMFITprofiles_190423_FIT.nc"

    rcp = pd.read_csv(rcp_path, delimiter="\t")

    # Load gfile to map to psin.
    with open(gfile_path, "rb") as f:
        gfile = pickle.load(f)
    R = gfile["R"]
    Z = gfile["Z"]
    Rs, Zs = np.meshgrid(R, Z)
    psin = gfile["PSIRZ_NORM"]

    rcp_coord = (rcp["R(cm)"].to_numpy() / 100 + rcp_shift, np.full(len(rcp["R(cm)"].values), -0.185))
    rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), rcp_coord)

    rcp_ne = rcp["Ne(E18 m-3)"].to_numpy() * 1e18

    # Load the OMFITprofiles data for getting ne at the separatrix.
    omfit_prof = netCDF4.Dataset(omfit_path)
    core_fit_psin = omfit_prof["psi_n"][:].data
    core_fit_ne = omfit_prof["n_e"][:].mean(axis=0).data
    core_fit_ne_err = np.sqrt(
        np.square(omfit_prof["n_e__uncertainty"][:].data).sum(axis=0))  # All zeros for some reason.
    f_core_fit_ne = interp1d(core_fit_psin, core_fit_ne)
    f_core_fit_ne_err = interp1d(core_fit_psin, core_fit_ne_err)
    nesep = f_core_fit_ne(1.0)

    rcp_r = rcp["R(cm)"].to_numpy() / 100 + rcp_shift
    return {"rcp_psin": rcp_psin, "rcp_ne": rcp_ne, "nesep": nesep, "rcp_r":rcp_r}


data195 = load_data(167195)
data423 = load_data(190423)
psin195 = data195["rcp_psin"][:-4]
psin423 = data423["rcp_psin"]
r195 = data195["rcp_r"][:-4]
r423 = data423["rcp_r"]
nenorm195 = data195["rcp_ne"][:-4] / data195["nesep"]
nenorm423 = data423["rcp_ne"] / data423["nesep"]


# Exponential fits to extract lambda_ne.
def exp_fit(x, a, b, c):
    return a * np.exp(-x / b) + c


mask195 = np.logical_and(psin195 > 1.0, psin195 < 1.23)
mask423 = np.logical_and(psin423 > 1.0, psin423 < 1.19)
popt195, pcov195 = curve_fit(exp_fit, r195[mask195]-2.3, nenorm195[mask195])
lamb195 = popt195[1]
popt423, pcov423 = curve_fit(exp_fit, r423[mask423]-2.3, nenorm423[mask423], maxfev=5000)
lamb423 = popt423[1]

c195 = "tab:purple"
c423 = "tab:red"
fig, ax1 = plt.subplots(figsize=(4, 3.5))

ax1.axvline(1.0, color="k", linestyle="--")
ax1.fill_between(psin195, nenorm195 * 0.8, nenorm195 * 1.2, alpha=0.3, color=c195)
ax1.fill_between(psin423, nenorm423 * 0.8, nenorm423 * 1.2, alpha=0.3, color=c423)
ax1.plot(psin195, nenorm195, label=167195, color=c195, lw=2)
ax1.plot(psin423, nenorm423, label=190423, color=c423, lw=2)
#ax1.plot(psin195[mask195], exp_fit(psin195[mask195], *popt195), color=c195, linestyle="--")
#ax1.plot(psin423[mask423], exp_fit(psin423[mask423], *popt423), color=c423, linestyle="--")
ax1.legend()
ax1.set_yscale("log")
ax1.set_ylim(3e-2, 2)
ax1.set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=12)
ax1.set_ylabel(r"$\mathdefault{n_e / n_e^{sep}}$", fontsize=12)
ax1.text(0.4, 0.6, r"$\mathdefault{\lambda_n} =$"+" {:.1f} cm".format(lamb195*100), transform=ax1.transAxes, color=c195, fontsize=12, rotation=-35)
ax1.text(0.12, 0.3, r"$\mathdefault{\lambda_n} =$"+" {:.1f} cm".format(lamb423*100), transform=ax1.transAxes, color=c423, fontsize=12, rotation=-50)

fig.tight_layout()
fig.show()
