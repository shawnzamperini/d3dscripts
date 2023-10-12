# Quickish script to calculate sigma = netarg / neomp using the LP and RCP values.
import get_lp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from scipy.interpolate import griddata, interp1d
import oedge_plots


# Load outer target LP data.
lpdict = get_lp.plot_lps(167195, 4000, 5000, "psin", bins=50, tunnel=False, showplot=False)
lp_psin = []
lp_ne = []
for i in range(0, len(lpdict["psin"])):
    pnum = int(lpdict["pnames"][i].split(" ")[1])
    if pnum >= 23 and pnum <= 53:
        lp_psin.append(lpdict["psin"][i])
        lp_ne.append(lpdict["ne (cm-3)"][i] * 1e6)
sort_idx = np.argsort(lp_psin)
lp_psin = np.array(lp_psin)[sort_idx]
lp_ne = np.array(lp_ne)[sort_idx]

# Load RCP data.
rcp_path = "/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP167195_2.tab"
rcp = pd.read_csv(rcp_path, delimiter="\t")
rcp_r = rcp["R(cm)"].to_numpy() / 100  # cm to m
rcp_ne = rcp["Ne(E18 m-3)"].to_numpy() * 1e18

# Load gfile to map RCP data to psin.
gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/167196/167196_3500.pickle"
with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)
R = gfile["R"]
Z = gfile["Z"]
Rs, Zs = np.meshgrid(R, Z)
psin = gfile["PSIRZ_NORM"]
rcp_coord = zip(rcp_r, np.full(len(rcp), -0.185))
rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))

# Can also include OEDGE results for comparison.
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-bg-shifted-ring-entry-09.nc"
op = oedge_plots.OedgePlots(ncpath)
div_rcp = op.fake_probe(2.215, 2.33, -0.185, -0.185, data="ne", plot="psin", show_plot=False)
div_psin = np.array(div_rcp["psin"])
div_ne = np.array(div_rcp["ne"])
ndsin = int(op.nc["NDSIN"][:].data) - 1  # First index of inner target (-1 for 0-indexing)
targ_ne = op.nc["KNDS"][:ndsin]
targ_ir = op.nc["IRDS"][:ndsin]  # rings of each target element.[:ndsin]
targ_psin = np.array([float(op.nc["PSIFL"][ir-1, 0].data) for ir in targ_ir])

# Sort them, make an interpolation function.
sort_idx = np.argsort(targ_psin)
targ_psin = targ_psin[sort_idx]
targ_ne = targ_ne[sort_idx]
targ_ir = targ_ir[sort_idx]
f_div_ne = interp1d(div_psin, div_ne)
f_targ_ne = interp1d(targ_psin, targ_ne)
minx = max(targ_psin.min(), div_psin.min())
maxx = min(targ_psin.max(), div_psin.max())
sigmax = np.linspace(minx, maxx, 100)
sigma = f_targ_ne(sigmax) / f_div_ne(sigmax)

# Just copied over from calc_sigma.py
psin_sigma_167195 = [1.19502619, 1.15942215, 1.12446229, 1.09153329]
sigma_167195 = [0.27281351,  1.09117605, 15.29851034,  5.50378404]

# Standard deviations for an error band.
err_max = 4 * sigma
err_min = 0.25 * sigma

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

ax1.axvline(1.0, color="k")
ax1.plot(lp_psin, lp_ne, label="LP")
ax1.plot(rcp_psin, rcp_ne, label="RCP")
ax1.plot(div_psin, div_ne, label="OEDGE - OMP")
ax1.plot(targ_psin, targ_ne, label="OEDGE - Target")
ax1.legend()
ax1.set_xlabel("Psin", fontsize=14)
ax1.set_ylabel("ne (m-3)", fontsize=14)
ax1.set_ylim(0, 0.4e20)

ax2.axhline(1.0, color="k", linestyle="--")
ax2.axvline(1.0, color="k", linestyle="--")
ax2.fill_between(sigmax, sigma / 2, sigma * 2, color="tab:red", alpha=0.3)
ax2.plot(sigmax, sigma, color="tab:red", lw=3)
ax2.scatter(psin_sigma_167195, sigma_167195, marker="^", s=75, color="tab:red", edgecolors="k", zorder=15)
ax2.set_xlabel(r"$\psi_n$", fontsize=14)
ax2.set_ylabel(r"$\sigma=\mathdefault{n_e^{targ}}$/$\mathdefault{n_e^{OMP}}$", fontsize=14)
ax2.set_ylim([0, None])

fig.tight_layout()
fig.show()

