import flan_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skewnorm
from tqdm import tqdm
from matplotlib.ticker import MultipleLocator, FuncFormatter


# Constants
amu_to_kg = 1.66e-27
elec = -1.602e-19

# Paths to each PT/NT simulation
pt_ncpath = '/pscratch/sd/z/zamp/flandir/pt_v5/sol_source_v1.nc'
pt_fp = flan_plots.FlanPlots(pt_ncpath)
nt_ncpath = '/pscratch/sd/z/zamp/flandir/nt_v5/sol_source_v1.nc'
nt_fp = flan_plots.FlanPlots(nt_ncpath)

# Frame to average over
tidx_min = -20
tidx_max = -1
pt_tmin = pt_fp.nc["geometry"]["time"][tidx_min]
pt_tmax = pt_fp.nc["geometry"]["time"][tidx_max]
nt_tmin = nt_fp.nc["geometry"]["time"][tidx_min]
nt_tmax = nt_fp.nc["geometry"]["time"][tidx_max]
print("Time ranges averaged over (should match!)")
print(f"  PT: {pt_tmin:.2e}  {pt_tmax:.2e}")
print(f"  NT: {nt_tmin:.2e}  {nt_tmax:.2e}")

# Radial index just OUTSIDE the LCFS. Important, the rest assumes this
# is the case
#xidx = 128  # For nt_v2
xidx = 80 # For nt_v3
xloc = pt_fp.nc["geometry"]["x"][xidx]
print(f"x = {xloc}")

# Radial coordinate and poloidal angle (same for both simulations)
r = pt_fp.nc["geometry"]["x"][:]
theta = pt_fp.nc["geometry"]["z"][:]

# Flag to average over the y coordinate or not
y_avg = True

# R, Z of the plasma center (doesn't need to be exact)
center_R = 1.8
center_Z = 0.0

def load_perp_drift(fp, drift, y_avg=False):

	# Midpoint of the y direction if y_avg = False
	yidx = int(len(fp.nc["geometry"]["y"][:]) / 2)


	# Grid coordinates
	grid_X = fp.nc["geometry"]["grid_X"][:]
	grid_Y = fp.nc["geometry"]["grid_Y"][:]
	grid_Z = fp.nc["geometry"]["grid_Z"][:]

	# Grad correct function for drift calculation
	if drift == "exb":
		drift_func = fp.calc_exb_drift
	elif drift == "gradb":
		drift_func = fp.calc_gradb_drift
	elif drift == "polarization":
		drift_func = fp.calc_polarization_drift
	elif drift == "curvature":
		drift_func = fp.calc_curvature_drift
	
	# Use "actual" for the actual velocity measured by Flan simulation
	elif drift == "actual":
		pass
	else:
		print("Error! Unknown option: {}".format(drift))

	# Average charge values and time average of those averages excluding
	# cells with zero counts.
	#qz = fp.nc["qz"][tidx_min:tidx_max].mean(axis=0)
	qz = fp.nc["output"]["qz"][tidx_min:tidx_max]
	Nz = fp.nc["output"]["Nz"][tidx_min:tidx_max]

	# Mask zeros
	mask = Nz != 0

	# Sum and count non-zero values along axis 0
	sum_ = np.sum(np.where(mask, qz, 0), axis=0)
	count = np.sum(mask, axis=0)

	# Compute mean safely
	qz = np.divide(sum_, count, where=count != 0)

	# Either y-average or index at a y location. Result is 1D array 
	# of len(theta)
	if y_avg:
		qz_avg = qz[xidx, :, :].mean(axis=1)
	else:
		qz_avg = qz[xidx, yidx, :]

	# Get drift components, assigning charge if necessary
	if drift == "exb":
		v_drift_X, v_drift_Y, v_drift_Z = drift_func()
	
	# Not actually drift, this is the simulation results
	elif drift == "actual":
		v_drift_X = fp.nc["output"]["v_X"][:]
		v_drift_Y = fp.nc["output"]["v_Y"][:]
		v_drift_Z = fp.nc["output"]["v_Z"][:]
	else:

		# Will need to load in z-loop so we can pass correct charge
		#v_drift_X, v_drift_Y, v_drift_Z = drift_func(charge)
		pass
	
	# Restrict to indicated x, y locations, or average over y. Can do this 
	# before the z-loop below for ExB and actual since they do not depend on
	# on charge
	if drift in ["exb", "actual"]:
		if y_avg:
			v_drift_X = v_drift_X[:, xidx, :, :].mean(axis=1)
			v_drift_Y = v_drift_Y[:, xidx, :, :].mean(axis=1)
			v_drift_Z = v_drift_Z[:, xidx, :, :].mean(axis=1)
		else:
			v_drift_X = v_drift_X[:, xidx, yidx, :]
			v_drift_Y = v_drift_Y[:, xidx, yidx, :]
			v_drift_Z = v_drift_Z[:, xidx, yidx, :]
		
		# Average in time for desired time window. Result is a 1D array with
		# length of z dimension
		v_avg_X = np.mean(v_drift_X[tidx_min:tidx_max], axis=0)
		v_avg_Y = np.mean(v_drift_Y[tidx_min:tidx_max], axis=0)
		v_avg_Z = np.mean(v_drift_Z[tidx_min:tidx_max], axis=0)

	# Since xidx is the first index outside the LCFS, the indices in grid
	# at xidx are those along the LCFS
	v_perp_X = np.zeros(len(theta))
	v_perp_Y = np.zeros(len(theta))
	v_perp_Z = np.zeros(len(theta))
	v_rad_R = np.zeros(len(theta))
	v_rad_Z = np.zeros(len(theta))
	direction_scalars = np.ones(len(theta))
	for zidx in tqdm(range(len(theta))):

		# Construct LCFS vector from its starting/ending coordinates. This
		# vector is going counter-clockwise in the poloidal direction
		lcfs_X0 = grid_X[xidx, yidx, zidx]
		lcfs_X1 = grid_X[xidx, yidx, zidx+1]
		lcfs_Y0 = grid_Y[xidx, yidx, zidx]
		lcfs_Y1 = grid_Y[xidx, yidx, zidx+1]
		lcfs_Z0 = grid_Z[xidx, yidx, zidx]
		lcfs_Z1 = grid_Z[xidx, yidx, zidx+1]
		lcfs_X = lcfs_X1 - lcfs_X0
		lcfs_Y = lcfs_Y1 - lcfs_Y0
		lcfs_Z = lcfs_Z1 - lcfs_Z0
		lcfs_R0 = np.sqrt(np.square(lcfs_X0) + np.square(lcfs_Y0))
		lcfs_R1 = np.sqrt(np.square(lcfs_X1) + np.square(lcfs_Y1))

		# Magnitude of LCFS segment (squared)
		lcfs_sq = np.square(lcfs_X) + np.square(lcfs_Y) + np.square(lcfs_Z)

		# If not doing ExB or actual values, need to load drifts here each 
		# loop iteration so we can pass the corresponding charge in.
		if drift not in ["exb", "actual"]:
			v_drift_X, v_drift_Y, v_drift_Z = drift_func(qz_avg[zidx])
			if y_avg:
				v_drift_X = v_drift_X[:, xidx, :, :].mean(axis=1)
				v_drift_Y = v_drift_Y[:, xidx, :, :].mean(axis=1)
				v_drift_Z = v_drift_Z[:, xidx, :, :].mean(axis=1)
			else:
				v_drift_X = v_drift_X[:, xidx, yidx, :]
				v_drift_Y = v_drift_Y[:, xidx, yidx, :]
				v_drift_Z = v_drift_Z[:, xidx, yidx, :]
		
			# Average in time for desired time window. Result is a 1D array with
			# length of z dimension
			v_avg_X = np.mean(v_drift_X[tidx_min:tidx_max], axis=0)
			v_avg_Y = np.mean(v_drift_Y[tidx_min:tidx_max], axis=0)
			v_avg_Z = np.mean(v_drift_Z[tidx_min:tidx_max], axis=0)

		# First calculate the parallel projection of the drift onto the LCFS
		scalar_proj = (v_avg_X[zidx] * lcfs_X + v_avg_Y[zidx] * lcfs_Y 
			+ v_avg_Z[zidx] * lcfs_Z)
		v_par_X = scalar_proj * lcfs_X / lcfs_sq
		v_par_Y = scalar_proj * lcfs_Y / lcfs_sq
		v_par_Z = scalar_proj * lcfs_Z / lcfs_sq

		# Then the perpendicular drift is the original vector minus the
		# parallel projection
		v_perp_X[zidx] = v_avg_X[zidx] - v_par_X
		v_perp_Y[zidx] = v_avg_Y[zidx] - v_par_Y
		v_perp_Z[zidx] = v_avg_Z[zidx] - v_par_Z
		v_perp_R = np.sqrt(np.square(v_perp_X[zidx])+np.square(v_perp_Y[zidx]))

		# Not done yet. We have the perpendicular drift, but we specifically
		# are interested in the portion perpendicular to the LCFS 
		# (alternatively, the portion parallel to the normal of the LCFS). So
		# first calculate the normal of the LCFS in (R,Z) by flipping the Z
		# coordinate. Do that at the midpoint of LCFS segment
		lcfs_mid_R = lcfs_R0 + 0.5 * (lcfs_R1 - lcfs_R0)
		lcfs_mid_Z = lcfs_Z0 + 0.5 * (lcfs_Z1 - lcfs_Z0)
		lcfs_norm_R = lcfs_mid_R
		lcfs_norm_Z = -lcfs_mid_Z
		
		# We want the normals to face radially outwards from the core center.
		# Just need a vector heading in that general direction.
		from_center_R = lcfs_mid_R - center_R
		from_center_Z = lcfs_mid_Z - center_Z

		# Dot product of normal with outward facing vector. If this is positive,
		# the normal is already facing outwards and vice versa.
		from_center_proj = from_center_R * lcfs_norm_R \
			+ from_center_Z * lcfs_norm_Z

		# Make sure normal faces outward now that we know which way. If it's
		# facing inwards, then just flip the components
		if (from_center_proj < 0): 
			lcfs_norm_R *= -1
			lcfs_norm_Z *= -1
		lcfs_norm_sq = np.square(lcfs_norm_R) + np.square(lcfs_norm_Z)

		# Now project the perpendicular velocity onto the outward normal
		# vector for the radially outward drift velocity in (R, Z)
		# First calculate the parallel projection of the drift onto the LCFS
		scalar_proj = v_perp_R * lcfs_norm_R \
			+ v_perp_Z[zidx] * lcfs_norm_Z
		v_rad_R[zidx] = scalar_proj * lcfs_norm_R / lcfs_norm_sq
		v_rad_Z[zidx] = scalar_proj * lcfs_norm_Z / lcfs_norm_sq

		# If scalar_proj is positive, then the vector is likewise facing
		# outward and direction_scalar = 1. Vice-versa.
		if (scalar_proj < 0): direction_scalars[zidx] = -1

	# Return magnitude of radial (outward) drift
	return np.sqrt(np.square(v_rad_R) + np.square(v_rad_Z)) * direction_scalars

# Actual perpendicular velocity
pt_v_perp_actual = load_perp_drift(pt_fp, "actual", y_avg)
nt_v_perp_actual = load_perp_drift(nt_fp, "actual", y_avg)

lw = 3
fontsize = 16
fig, ax1 = plt.subplots(1, 1, figsize=(5, 4))

for ax in [ax1]:
	ax.tick_params(axis="both", which="both", labelsize=fontsize-2)
	ax.set_xlim(-np.pi, np.pi)

	# Set tick locations at multiples of π/2
	ax.xaxis.set_major_locator(MultipleLocator(base=np.pi / 2))

	# Custom formatter to show π symbols
	def format_pi(x, pos):
		labels = {
			-np.pi: r'$\mathdefault{-\pi}$',
			-np.pi/2: r'$\mathdefault{-\pi / 2}$',
			0: r'$0$',
			np.pi/2: r'$\mathdefault{\pi / 2}$',
			np.pi: r'$\mathdefault{\pi}$'
		}
		return labels.get(x, f'{x:.2f}')

	ax.xaxis.set_major_formatter(FuncFormatter(format_pi))

ax1.plot(theta, pt_v_perp_actual, label="PT")
ax1.plot(theta, nt_v_perp_actual, label="NT")
ax1.set_ylabel("Speed Leaving Core (m/s)", fontsize=fontsize) 
ax1.legend()
fig.tight_layout()
fig.show()
