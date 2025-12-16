import pickle
import numpy as np
import matplotlib.pyplot as plt
from S3X_EIRENE import S3X_EIRENE


# Load pickled TS and RCP data from 167196 and the associated diagnostic shots
with open("/global/homes/z/zamp/data/exp_data_167196.pickle", "rb") as f:
	exp_data = pickle.load(f)

# Load an S3XE case
simulation_folder = "/pscratch/sd/z/zamp/soledge3x_rundir/diiid_167196_v1"

sim = S3X_EIRENE()
sim.load_S3XE_simulation(
    base_folder             = simulation_folder,
    run_dir                 = "run_dir", # Default
    plasma_folder           = "Plasma", # Default
    load_restart            = False, # Default
    load_final              = False, # Default
    load_only_indices       = None, # Default, all plasma files. can be integer, or list/range (eg. range(0,1000,100)), put -1 to load only last file
    load_emergency_saves    = False, # Default
    force_reload            = False
)

# Load QuadPlasma object
qp = sim.SOLEDGE_plasma

# Load QuadFields
n_e = sim.SOLEDGE_plasma.species["e-"]["n"]
T_e = sim.SOLEDGE_plasma.species["e-"]["T"]
T_i = sim.SOLEDGE_plasma.species["D+"]["T"]

# Coordinates along the TS core system line
ts_core_line_npoints = 100
ts_core_r = exp_data["ts_core"]["r (m)"][0]
ts_core_line_r = np.full(ts_core_line_npoints, ts_core_r)
ts_core_line_z = np.linspace(0.53, 1.0, ts_core_line_npoints)

# Density and temperature along TS core system for final plasma frame
s3xe_ne_ts_core = qp.get_values_from_RZ(n_e, 
	list(zip(ts_core_line_r, ts_core_line_z)), -1)
s3xe_te_ts_core = qp.get_values_from_RZ(T_e, 
	list(zip(ts_core_line_r, ts_core_line_z)), -1)

# Similarly for along the RCP location
rcp_line_npoints = 100
rcp_line_r = np.linspace(2.25, 2.37, rcp_line_npoints)
rcp_line_z = np.full(rcp_line_npoints, -0.188) # RCP at Z = -0.188 m
s3xe_ne_rcp = qp.get_values_from_RZ(n_e, 
	list(zip(rcp_line_r, rcp_line_z)), -1)
s3xe_te_rcp = qp.get_values_from_RZ(T_e, 
	list(zip(rcp_line_r, rcp_line_z)), -1)

# Flatten TS data so we can plot all the data points
ts_core_ntimes = exp_data["ts_core"]["ne (m-3)"].shape[1]
ts_core_z_flat = np.repeat(exp_data["ts_core"]["z (m)"], ts_core_ntimes)
ts_core_ne_flat = exp_data["ts_core"]["ne (m-3)"].flatten()
ts_core_te_flat = exp_data["ts_core"]["te (eV)"].flatten()

# RCP data, put into a single array
rcp_r = np.array([])
rcp_z = np.array([])
rcp_ne = np.array([])
rcp_te = np.array([])
rcp_m = np.array([])
rcp_vp = np.array([])
for shot in range(167192, 167196):
	for plunge in [1, 2]:
		rcp = exp_data["rcp"]["MP{}_{}".format(shot, plunge)]
		rcp_r = np.append(rcp_r, np.array(rcp["R(cm)"]) / 100)
		rcp_z = np.append(rcp_z, rcp["Z(m)"])
		rcp_ne = np.append(rcp_ne, np.array(rcp["Ne(E18 m-3)"]) * 1e18)
		rcp_te = np.append(rcp_te, rcp["Te(eV)"])
		rcp_m = np.append(rcp_m, rcp["Machn"])
		rcp_vp = np.append(rcp_vp, rcp["Vp (V)"])

# Separatrix location along each diagnostic
ts_core_sep_z = 0.733
rcp_sep_r = 2.21667

fontsize = 16
fig, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2, 3, figsize=(15, 8))

# Density at TS core location
ax1.axvline(ts_core_sep_z, color="k", linestyle="--")
ax1.scatter(ts_core_z_flat, ts_core_ne_flat, label="TS", alpha=0.3, 
	c="tab:green", edgecolors="k", s=50)
ax1.plot(ts_core_line_z, s3xe_ne_ts_core, lw=3, color="k")
ax1.plot(ts_core_line_z, s3xe_ne_ts_core, lw=2, color="tab:green", label="S3XE")
ax1.legend()
ax1.set_xlabel("Z (m)", fontsize=fontsize)
ax1.set_ylabel("ne (m-3)", fontsize=fontsize)
ax1.set_xlim([0.4, 1.1])
ax1.set_ylim([0.0, 4e19])

# Temperature at TS core location
ax3.axvline(ts_core_sep_z, color="k", linestyle="--")
ax3.scatter(ts_core_z_flat, ts_core_te_flat, label="TS", alpha=0.3, 
	c="tab:red", edgecolors="k", s=50)
ax3.plot(ts_core_line_z, s3xe_te_ts_core, lw=3, color="k")
ax3.plot(ts_core_line_z, s3xe_te_ts_core, lw=2, color="tab:red", label="S3XE")
ax3.legend()
ax3.set_xlabel("Z (m)", fontsize=fontsize)
ax3.set_ylabel("Te (eV)", fontsize=fontsize)
ax3.set_xlim([0.4, 1.1])
ax3.set_ylim([0.0, 1100])

# Density at RCP location
ax2.axvline(rcp_sep_r, color="k", linestyle="--")
ax2.scatter(rcp_r, rcp_ne, label="RCP", alpha=0.3, c="tab:green", 
	edgecolors="k", s=50)
ax2.plot(rcp_line_r, s3xe_ne_rcp, lw=3, color="k")
ax2.plot(rcp_line_r, s3xe_ne_rcp, lw=2, color="tab:green")
ax2.legend()
ax2.set_xlabel("R (m)", fontsize=fontsize)
ax2.set_ylabel("ne (m-3)", fontsize=fontsize)
ax2.set_ylim([0.0, 2.5e19])

# Temperature at RCP location
ax4.axvline(rcp_sep_r, color="k", linestyle="--")
ax4.scatter(rcp_r, rcp_te, label="RCP", alpha=0.3, c="tab:red", 
	edgecolors="k", s=50)
ax4.plot(rcp_line_r, s3xe_te_rcp, lw=3, color="k")
ax4.plot(rcp_line_r, s3xe_te_rcp, lw=2, color="tab:red")
ax4.legend()
ax4.set_xlabel("R (m)", fontsize=fontsize)
ax4.set_ylabel("Te (eV)", fontsize=fontsize)
ax4.set_ylim([0.0, 90])

# Mach number at RCP location
ax5.axvline(rcp_sep_r, color="k", linestyle="--")
ax5.axhline(0.0, color="k")
ax5.scatter(rcp_r, rcp_m, label="RCP", alpha=0.3, c="tab:purple", 
	edgecolors="k", s=50)
#ax5.plot(rcp_line_r, s3xe_te_rcp, lw=3, color="k")
#ax5.plot(rcp_line_r, s3xe_te_rcp, lw=2, color="tab:red")
ax5.legend()
ax5.set_xlabel("R (m)", fontsize=fontsize)
ax5.set_ylabel("Mach", fontsize=fontsize)
ax5.set_ylim([-1, 1])

# Plasma potential number at RCP location
ax6.axvline(rcp_sep_r, color="k", linestyle="--")
ax6.scatter(rcp_r, rcp_vp, label="RCP", alpha=0.3, c="tab:orange", 
	edgecolors="k", s=50)
#ax6.plot(rcp_line_r, s3xe_te_rcp, lw=3, color="k")
#ax6.plot(rcp_line_r, s3xe_te_rcp, lw=2, color="tab:red")
ax6.legend()
ax6.set_xlabel("R (m)", fontsize=fontsize)
ax6.set_ylabel("Vp (V)", fontsize=fontsize)
ax6.set_ylim([0.0, None])

fig.tight_layout()
fig.show()
