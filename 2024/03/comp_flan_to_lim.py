# Goal here is to make a plot comparing the inferred W density profile 
# from 3DLIM to that from Flan.
import flan_plots
import LimPlots
import numpy as np
import matplotlib.pyplot as plt


# COpying this block form gkyl_imp_sim for 3DLIm data.
lim_path = "/Users/zamperini/Documents/d3d_work/lim_runs/167196/167196-a2-tor240-blob-013-018d-noprobe.nc"
lp = LimPlots.LimPlots(lim_path)
lpdata = lp.plot_par_rad("nz", 21, charge="all")
mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 15
rorigin = 2.282
rad_locs = rorigin - lpdata["Y"][0]
absfac = 1e15
lim_nz = lpdata["Z"][mididx].data * absfac
lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))

# Load standard plot of the Flan run.
flan_path = "/Users/zamperini/flandir/reg_testcase1/saved_results_4/reg_testcase1.nc"
fp = flan_plots.FlanPlots(flan_path)
fp_data = fp.plot_profiles(["ne", "nz"], plot_z=0.3125, normtype=["log","log"], vmin=[1e18, 1], vmax=[2e19, 1e3], skip_video=True)
x = fp_data["x"]

# Impurity density averaged over the y coordinate.
nz = fp_data["data"][1].mean(axis=2)

# Pick a range of frames after some time to equlibriate.
nz_avg = nz[75:].mean(axis=0)

# Normalize to the value at R=2.33, a value that both simulations overlap.
match_r = 2.34
#match_rmrs = 0.08
#rcp_rsep = 2.216
# print("match_rmrsep = {:.3f}".format(match_rmrs - gk_rsep))
#match_lim_idx = np.argmin(np.abs(rad_locs-rcp_rsep - match_rmrs))
match_lim_idx = np.argmin(np.abs(rad_locs - match_r))
lim_nz_norm = lim_nz / lim_nz[match_lim_idx]
#match_gkyl_idx = np.argmin(np.abs(plasma.r-gk_rsep - match_rmrs))
#rad_imp_dens_norm = rad_imp_dens / rad_imp_dens[match_gkyl_idx]
match_flan_idx = np.argmin(np.abs(x - match_r))
nz_avg_norm = nz_avg / nz_avg[match_flan_idx]

fig, ax1 = plt.subplots()
ax1.plot(rad_locs, lim_nz_norm)
ax1.plot(x, nz_avg_norm)
ax1.legend()
ax1.set_yscale("log")
fig.tight_layout()
fig.show()
