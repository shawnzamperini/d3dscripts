import pickle
import matplotlib.pyplot as plt
import numpy as np


# Load pickled data for #167196. This includes data from #167192-167195 since
# these were repeat diagnostic shots. This means:
#	- The strike point was swept to fill out LP profiles
#   - RCP on MiMES to get midplane ne, Te, Mach, etc. data. A collector probe
#       was on MiMES for 167196.
with open("/home/zamp/data/exp_data_167196.pickle", "rb") as f:
	exp = pickle.load(f)

# Plotting 
ts_color = "tab:red"
lp_color = "tab:purple"
rcp_color = "tab:cyan"
fontsize = 14
fig, axs = plt.subplots(3, 3, figsize=(9, 8))

# -----------------------------------
#       Thomson scattering plots. 
# -----------------------------------

# Get x-axis data first. I mask the data to a time 
# range I have identified as during the density flat top.
t = exp["ts_core"]["t (ms)"]
tmin = exp["ts_core"]["t_min (ms)"]
tmax = exp["ts_core"]["t_max (ms)"]
mask = np.logical_and(t >= tmin, t <= tmax)

# Loop through each TS chord. Chords are vertically stacked along Z. 
axs[0][0].axvline(0.73, color="k", linestyle="--")
axs[0][1].axvline(0.73, color="k", linestyle="--")
for chord in range(exp["ts_core"]["te (eV)"].shape[0] - 1):

	# Extract signals, masking to flat top range
	z = np.repeat(exp["ts_core"]["z (m)"][chord], sum(mask))
	te = exp["ts_core"]["te (eV)"][chord][mask]
	ne = exp["ts_core"]["ne (m-3)"][chord][mask]

	# Te
	axs[0][0].set_title("Thomson Scattering", color=ts_color)
	axs[0][0].scatter(z, te, color=ts_color, alpha=0.2, s=25)
	axs[0][0].set_xlabel("Z (m)", fontsize=fontsize)
	axs[0][0].set_ylabel("Te (eV)", fontsize=fontsize)

	# ne
	axs[0][1].set_title("Thomson Scattering", color=ts_color)
	axs[0][1].scatter(z, ne, color=ts_color, alpha=0.2, s=25)
	axs[0][1].set_xlabel("Z (m)", fontsize=fontsize)
	axs[0][1].set_ylabel("ne (m-3)", fontsize=fontsize)

# ----------------------------------------------
#       Reciprocating Langmuir probe plots. 
#
# RCP is located at Z = -0.188 m
# ----------------------------------------------

# Loop through each plunge, 2 per shot for #167192-5
for plunge in exp["rcp"].keys():
	
	r = np.array(exp["rcp"][plunge]["R(cm)"]) / 100  # cm to m
	te = np.array(exp["rcp"][plunge]["Te(eV)"])
	ne = np.array(exp["rcp"][plunge]["Ne(E18 m-3)"]) * 1e18  # to m-3
	mach = np.array(exp["rcp"][plunge]["Machn"])
	qpar = np.array(exp["rcp"][plunge]["Q_par(MW/m^2)"])

	# Spikes in data at the most inward location are due to probe tip arcing
	axs[0][2].scatter(r, te, color=rcp_color)
	axs[0][2].set_xlabel("R (m)", fontsize=fontsize)
	axs[0][2].set_ylabel("Te (eV)", fontsize=fontsize)
	axs[0][2].set_title("RCP", fontsize=fontsize, color=rcp_color)

	# Spikes in data at the most inward location are due to probe tip arcing
	axs[1][0].scatter(r, ne, color=rcp_color)
	axs[1][0].set_xlabel("R (m)", fontsize=fontsize)
	axs[1][0].set_ylabel("ne (m-3)", fontsize=fontsize)
	axs[1][0].set_ylim([-0.5e18, 2e19])
	axs[1][0].set_title("RCP", fontsize=fontsize, color=rcp_color)

	axs[1][1].scatter(r, mach, color=rcp_color)
	axs[1][1].set_xlabel("R (m)", fontsize=fontsize)
	axs[1][1].set_ylabel("Mach Number", fontsize=fontsize)
	axs[1][1].set_ylim([-1, 0.6])
	axs[1][1].set_title("RCP", fontsize=fontsize, color=rcp_color)

	axs[1][2].scatter(r, qpar, color=rcp_color)
	axs[1][2].set_xlabel("R (m)", fontsize=fontsize)
	axs[1][2].set_ylabel(r"$\mathdefault{q_{||}\ (MW/m^2)}$", fontsize=fontsize)
	axs[1][2].set_ylim([0, 10])
	axs[1][2].set_title("RCP", fontsize=fontsize, color=rcp_color)

# ----------------------------------------------
#       Target Langmuir probe plots. 
# ----------------------------------------------

# Shelf (outer target) data. Can plot inner data if you care but there isn't
# much.
axs[2][0].scatter(exp["lp"]["shelf"]["psin"], exp["lp"]["shelf"]["Te (eV)"], 
	alpha=0.3, color=lp_color)
#axs[2][0].scatter(exp["lp"]["inner"]["psin"], exp["lp"]["inner"]["Te (eV)"], 
#	alpha=0.3, color="r")
axs[2][0].set_xlim([0.97, 1.11])
axs[2][0].set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=fontsize)
axs[2][0].set_ylabel("Te (eV)", fontsize=fontsize)
axs[2][0].set_title("Target LPs", fontsize=fontsize, color=lp_color)

axs[2][1].scatter(exp["lp"]["shelf"]["psin"], exp["lp"]["shelf"]["ne (m-3)"], 
	alpha=0.3, color=lp_color)
#axs[2][1].scatter(exp["lp"]["inner"]["psin"], exp["lp"]["inner"]["ne (m-3)"], 
#	alpha=0.3, color="r")
axs[2][1].set_xlim([0.97, 1.11])
axs[2][1].set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=fontsize)
axs[2][1].set_ylabel("ne (m-3)", fontsize=fontsize)
axs[2][1].set_title("Target LPs", fontsize=fontsize, color=lp_color)

axs[2][2].scatter(exp["lp"]["shelf"]["psin"], exp["lp"]["shelf"]["qperp (W/m2)"] / 1e6, 
	alpha=0.3, color=lp_color)
#axs[2][2].scatter(exp["lp"]["inner"]["psin"], exp["lp"]["inner"]["qperp (W/m2)"] / 1e6, 
#	alpha=0.3, color="r")
axs[2][2].set_xlim([0.97, 1.11])
axs[2][2].set_xlabel(r"$\mathdefault{\psi_n}$", fontsize=fontsize)
axs[2][2].set_ylabel(r"$\mathdefault{q}_{\perp}\ \mathdefault{(MW/m^2)}$", fontsize=fontsize)
axs[2][2].set_title("Target LPs", fontsize=fontsize, color=lp_color)

fig.tight_layout()
fig.show()
