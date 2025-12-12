import flan_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skewnorm
from tqdm import tqdm
from matplotlib.ticker import MultipleLocator, FuncFormatter


# Constants
amu_to_kg = 1.66e-27
elec = -1.602e-19

# Flan simulation. Below we assume that we can get a full poloidal profile
# by picking an arbitrary y value (say, the halfway point), pick a radial (x)
# coordiante, say at the separatrix, and then plot it against the poloidal
# angle (z, 0-2pi). 
#ncpath = "/home/zamp/flandir/nt_v2/nt_v2_coll.nc"
ncpath = "/home/zamp/flandir/nt_v3/nt_v3_coll.nc"
fp = flan_plots.FlanPlots(ncpath)

# Frame to average over
tidx_min = -30
tidx_max = -1

# Radial index just OUTSIDE the LCFS. Important, the rest assumes this
# is the case
#xidx = 128  # For nt_v2
xidx = 80 # For nt_v3

# Flag to average over the y coordinate or not
y_avg = True

# Midpoint of the y direction if y_avg = False
yidx = int(len(fp.nc["geometry"]["y"][:]) / 2)

# Radial coordinate and poloidal angle
r = fp.nc["geometry"]["x"][:]
theta = fp.nc["geometry"]["z"][:]

# Grid coordinates
grid_X = fp.nc["geometry"]["grid_X"][:]
grid_Y = fp.nc["geometry"]["grid_Y"][:]
grid_Z = fp.nc["geometry"]["grid_Z"][:]

# R, Z of the plasma center (doesn't need to be exact)
center_R = 1.8
center_Z = 0.0

def calc_char_time_E():
	"""
	Calculate the characteristic time scale for EX, EY and EZ, returning
	"""

	# Pull out some variables, take y average if asked
	time = fp.nc["geometry"]["time"][:]
	E_X = fp.nc["background"]["E_X"][:]
	E_Y = fp.nc["background"]["E_Y"][:]
	E_Z = fp.nc["background"]["E_Z"][:]
	E = np.sqrt(np.square(E_X) + np.square(E_Y) + np.square(E_Z))
	
	# Estimate time derivative as change from previous frame to current
	char_t_X = np.zeros(E_X.shape)
	char_t_Y = np.zeros(E_X.shape)
	char_t_Z = np.zeros(E_X.shape)
	char_t = np.zeros(E_X.shape)
	for t in range(1, E_X.shape[0]):

		# Calculate time derivative
		dt = time[t] - time[t-1]
		dEX_dt = (E_X[t] - E_X[t-1]) / dt
		dEY_dt = (E_Y[t] - E_Y[t-1]) / dt
		dEZ_dt = (E_Z[t] - E_Z[t-1]) / dt
		dE_dt = (E[t] - E[t-1]) / dt

		# Then characteristic time is E / dE/dt
		char_t_X[t] = np.abs(E_X[t] / dEX_dt)
		char_t_Y[t] = np.abs(E_Y[t] / dEY_dt)
		char_t_Z[t] = np.abs(E_Z[t] / dEZ_dt)
		char_t[t] = np.abs(E[t] / dE_dt)
	
	return char_t_X, char_t_Y, char_t_Z, char_t

def calc_char_length_E_perp():

	# Pull out components
	EX = fp.nc["background"]["E_X"][:]
	EY = fp.nc["background"]["E_Y"][:]
	EZ = fp.nc["background"]["E_Z"][:]
	BX = fp.nc["background"]["B_X"][:]
	BY = fp.nc["background"]["B_Y"][:]
	BZ = fp.nc["background"]["B_Z"][:]
	gradEX = fp.nc["background"]["gradE_X"][:]
	gradEY = fp.nc["background"]["gradE_Y"][:]
	gradEZ = fp.nc["background"]["gradE_Z"][:]

	# Calculate magnitudes
	Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

	# Calculate perpendicular to B components for electric field...
	scalar_proj = (EX * BX + EY * BY + EZ * BZ)
	EparX = scalar_proj * BX / Bsq
	EparY = scalar_proj * BY / Bsq
	EparZ = scalar_proj * BZ / Bsq
	EperpX = EX - EparX
	EperpY = EY - EparY
	EperpZ = EZ - EparZ

	# ...and again for the gradient
	scalar_proj = (gradEX * BX + gradEY * BY + gradEZ * BZ)
	gradEparX = scalar_proj * BX / Bsq
	gradEparY = scalar_proj * BY / Bsq
	gradEparZ = scalar_proj * BZ / Bsq
	gradEperpX = gradEX - gradEparX
	gradEperpY = gradEY - gradEparY
	gradEperpZ = gradEZ - gradEparZ

	# Magnitude of perpendicular components
	Eperp = np.sqrt(np.square(EperpX) + np.square(EperpY)
		+ np.square(EperpZ))
	gradEperp = np.sqrt(np.square(gradEperpX) + np.square(gradEperpY) 
		+ np.square(gradEperpZ))

	# Characteristic length is then just the absolute of the field divided
	# by it's gradient.
	char_l_Eperp = np.abs(Eperp / gradEperp)
	return char_l_Eperp
	

def calc_cycl_freq():

	# Pull out some arrays
	vX = fp.nc["output"]["v_X"][:]
	vY = fp.nc["output"]["v_Y"][:]
	vZ = fp.nc["output"]["v_Z"][:]
	BX = fp.nc["background"]["B_X"][:]
	BY = fp.nc["background"]["B_Y"][:]
	BZ = fp.nc["background"]["B_Z"][:]
	mz_kg = fp.nc["input"]["imp_mass_amu"][:] * amu_to_kg
	qz = fp.nc["output"]["qz"][:]

	# Magnetic field magnitude squared
	Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

	# Magnetic field unit vector
	bX = BX / Bsq
	bY = BY / Bsq
	bZ = BZ / Bsq
	
	# Cyclotron frequency
	omega_c = np.abs(qz * elec) * np.sqrt(Bsq) / mz_kg
	return omega_c

def calc_vperp():

	# Pull out some arrays
	vX = fp.nc["output"]["v_X"][:]
	vY = fp.nc["output"]["v_Y"][:]
	vZ = fp.nc["output"]["v_Z"][:]
	BX = fp.nc["background"]["B_X"][:]
	BY = fp.nc["background"]["B_Y"][:]
	BZ = fp.nc["background"]["B_Z"][:]

	# Magnetic field magnitude squared
	Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

	# Perpendicular to B velocity. Can obtain the perpendicular components
	# by subtracting the parallel projection of the impurity velocity from 
	# the total velocity vector.
	scalar_proj = (vX * BX + vY * BY + vZ * BZ)
	vparX = scalar_proj * BX / Bsq
	vparY = scalar_proj * BY / Bsq
	vparZ = scalar_proj * BZ / Bsq
	vperpX = vX - vparX
	vperpY = vY - vparY
	vperpZ = vZ - vparZ
	vperp = np.sqrt(np.square(vperpX) + np.square(vperpY) + np.square(vperpZ))
	
	return vperp

def load_epol(y_avg=False):
	
	E_X = fp.nc["background"]["E_X"][:]
	E_Y = fp.nc["background"]["E_Y"][:]
	E_Z = fp.nc["background"]["E_Z"][:]
	if y_avg:
		E_X = E_X[:, xidx, :, :].mean(axis=1)
		E_Y = E_Y[:, xidx, :, :].mean(axis=1)
		E_Z = E_Z[:, xidx, :, :].mean(axis=1)
	else:
		E_X = E_X[:, xidx, yidx, :]
		E_Y = E_Y[:, xidx, yidx, :]
		E_Z = E_Z[:, xidx, yidx, :]
	
	# Average in time for desired time window. Result is a 1D array with
	# length of z dimension
	E_avg_X = np.mean(E_X[tidx_min:tidx_max], axis=0)
	E_avg_Y = np.mean(E_Y[tidx_min:tidx_max], axis=0)
	E_avg_Z = np.mean(E_Z[tidx_min:tidx_max], axis=0)

	# Since xidx is the first index outside the LCFS, the indices in grid
	# at xidx are those along the LCFS
	E_pol_R = np.zeros(len(theta))
	E_pol_Z = np.zeros(len(theta))
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
		lcfs_R = lcfs_R1 - lcfs_R0

		# Magnitude of LCFS segment (squared)
		lcfs_sq = np.square(lcfs_X) + np.square(lcfs_Y) + np.square(lcfs_Z)

		# First calculate the parallel projection of E onto the LCFS
		scalar_proj = (E_avg_X[zidx] * lcfs_X + E_avg_Y[zidx] * lcfs_Y 
			+ E_avg_Z[zidx] * lcfs_Z)
		E_par_X = scalar_proj * lcfs_X / lcfs_sq
		E_par_Y = scalar_proj * lcfs_Y / lcfs_sq
		E_par_Z = scalar_proj * lcfs_Z / lcfs_sq
		E_par_R = np.sqrt(np.square(E_par_X) + np.square(E_par_Y))

		# We have the parallel E component, but we are after the poloidal
		# component, which is that that is parallel to the LCFS segment in
		# (R, Z).
		scalar_proj = E_par_R * lcfs_R + E_par_Z * lcfs_Z
		E_pol_R[zidx] = scalar_proj * lcfs_R / lcfs_sq
		E_pol_Z[zidx] = scalar_proj * lcfs_Z / lcfs_sq

		# Save the direction scalar. The sign ultimately depends on what 
		# direction we are traversing the LCFS, but it's not too important
		# which is which for our purposes at the moment.
		if (scalar_proj < 0): direction_scalars[zidx] = -1
	
	# Now return magnitude
	return np.sqrt(np.square(E_pol_R) + np.square(E_pol_Z)) * direction_scalars


def load_perp_drift(drift, y_avg=False):

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

		"""
		# We're not done yet. Now we need to determine if the vector is
		# pointing inwards or outwards from the core. Start by creating the
		# vector that points from the LCFS towards the center of the core.
		# We can do this in 2D (R,Z).
		lcfs_R0 = np.sqrt(np.square(lcfs_X0) + np.square(lcfs_Z0))
		lcfs_to_center_R = lcfs_R0 - center_R
		lcfs_to_center_Z = lcfs_Z0 - center_Z

		# Dot product between the velocity vector and vector pointing inwards
		inward_proj = v_perp_R * lcfs_to_center_R \
			+ v_perp_Z[zidx] * lcfs_to_center_Z

		# If inward proj is positive, then the velocity points inwards, and
		# vice-versa. Set direction_scalar = 1 for outward transport, and 
		# -1 for inward transport. Feels like that's the expected convention.
		direction_scalars[zidx] = 1
		if (inward_proj > 0): direction_scalars[zidx] = -1 
		"""
	
	"""
	# Then return the magnitude of it times the direction scalar (positive
	# = outward). 
	return np.sqrt(np.square(v_perp_X) + np.square(v_perp_Y) 
		+ np.square(v_perp_Z)) * direction_scalars
	"""
	# Return magnitude of radial (outward) drift
	return np.sqrt(np.square(v_rad_R) + np.square(v_rad_Z)) * direction_scalars


# Load each radial drift component (+ = outwards across LCFS)
print("ExB...")
v_perp_exb = load_perp_drift("exb", y_avg)
print("Grad-B...")
v_perp_gradb = load_perp_drift("gradb", y_avg)
print("Polarization...")
v_perp_pol = load_perp_drift("polarization", y_avg)
print("Curvature...")
v_perp_curv = load_perp_drift("curvature", y_avg)

# Total drift
v_drift_total = v_perp_exb + v_perp_gradb + v_perp_pol + v_perp_curv

# Actual perpendicular velocity
v_perp_actual = load_perp_drift("actual", y_avg)

# Poloidal electric field magnitude
Epol = load_epol(y_avg)

# Characteristic time for electric field fluctuations
char_t_X, char_t_Y, char_t_Z, char_t = calc_char_time_E()
char_l_Eperp = calc_char_length_E_perp()

# Characteristic times at the LCFS, average in time
if y_avg:
	char_t_X_lcfs = char_t_X[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
	char_t_Y_lcfs = char_t_Y[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
	char_t_Z_lcfs = char_t_Z[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
	char_t_lcfs = char_t[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
	char_l_Eperp_lcfs = char_l_Eperp[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
else:
	char_t_X_lcfs = char_t_X[tidx_min:tidx_max, xidx, yidx].mean(axis=0)
	char_t_Y_lcfs = char_t_Y[tidx_min:tidx_max, xidx, yidx].mean(axis=0)
	char_t_Z_lcfs = char_t_Z[tidx_min:tidx_max, xidx, yidx].mean(axis=0)
	char_t_lcfs = char_t[tidx_min:tidx_max, xidx, yidx].mean(axis=0)
	char_l_Eperp_lcfs = char_l_Eperp[tidx_min:tidx_max, xidx, yidx].mean(axis=0)

# Cyclotron frequency for impurities
omega_c = calc_cycl_freq()

# Cyclotron frequency at the LCFS
if y_avg:
	omega_c_lcfs = omega_c[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
else:
	omega_c_lcfs = omega_c[tidx_min:tidx_max, xidx, yidx].mean(axis=0)

# Average vperp for impurities (perpendicular to B, getting hairy here)
vperp_to_B = calc_vperp()

# Vperp to B at the LCFS
if y_avg:
	vperp_to_B_lcfs = vperp_to_B[tidx_min:tidx_max, xidx].mean(axis=0).mean(axis=0)
else:
	vperp_to_B_lcfs = vperp_to_B[tidx_min:tidx_max, xidx, yidx].mean(axis=0)

# Larmor radius at LCFS
larmor_radius_lcfs = vperp_to_B_lcfs / omega_c_lcfs

# Flattened subsets for plotting
char_t_X_sub = char_t_X[tidx_min:].flatten()
char_t_Y_sub = char_t_Y[tidx_min:].flatten()
char_t_Z_sub = char_t_Z[tidx_min:].flatten()
omega_c_sub = omega_c[tidx_min:].flatten()

# When this is >~ 0.1 the guiding-center approximation is not valid
criteria_X = char_t_X_sub * omega_c_sub
criteria_Y = char_t_Y_sub * omega_c_sub
criteria_Z = char_t_Z_sub * omega_c_sub
bins = np.logspace(np.log10(0.001), np.log10(1000), num=20)

# Citeria at the LCFS for total E magnitude, broken into "good" and "bad" regions
criteria_lcfs = char_t_lcfs * omega_c_lcfs
criteria_length_lcfs = char_l_Eperp_lcfs / larmor_radius_lcfs
good_mask = criteria_lcfs > 10
criteria_lcfs_good = np.copy(criteria_lcfs)
criteria_lcfs_good[good_mask] = np.nan
criteria_lcfs_bad = np.copy(criteria_lcfs)
criteria_lcfs_bad[~good_mask] = np.nan

# List of colors for good/bad indicators on the bars. Just manually assign
start_bad = 14
bar_colors = ["tab:green"] * (len(bins) - start_bad) + ["tab:red"] * start_bad

# Difference between Flan and the drifts
diff = v_perp_actual - v_drift_total

# Helper function to caluclate fraction of histogram above a value
def frac_above(limit, n, bins):
	mask = bins[:-1] > limit
	return n[mask].sum()

lw = 3
fontsize = 16
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13, 10))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Red gradient in background of criteria plot to give a sense of when you're
# getting into trouble
gradient = np.geomspace(1, 1e-5, 1000).reshape(-1, 1)
ax2.imshow(gradient, extent=[-np.pi, np.pi, 0.1, 20], 
	origin="lower", aspect="auto", cmap="Reds", alpha=0.5)

#for ax in [ax1, ax2, ax3, ax4]:
for ax in [ax1, ax2]:
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

ax1.axhline(0.0, color="k")
#ax1.plot(theta, v_perp_exb, label="ExB", lw=lw)
#ax1.plot(theta, v_perp_gradb, label="Grad-B", lw=lw)
#ax1.plot(theta, v_perp_pol, label="Polarization", lw=lw)
#ax1.plot(theta, v_perp_curv, label="Curvature", lw=lw)
ax1.plot(theta, v_drift_total, lw=lw+1, color="k")
ax1.plot(theta, v_drift_total, label="Total Drift", lw=lw, color="tab:red")
ax1.plot(theta, v_perp_actual, lw=lw+1, color="k")
ax1.plot(theta, v_perp_actual, label="Flan", lw=lw, color="tab:purple")
ax1.set_ylim([-15000, None])
ax1.legend(fontsize=fontsize-6, ncols=2, loc="lower center")
ax1.set_xlabel("Poloidal Angle", fontsize=fontsize)
ax1.set_ylabel("Speed Leaving Core (m/s)", fontsize=fontsize) 

ax2.axhline(10.0, color="k", linestyle="--")
ax2.plot(theta, criteria_lcfs, lw=lw+1, color="k")
ax2.plot(theta, criteria_lcfs, lw=lw, color="tab:pink", label=r"$\mathdefault{\tau_E \omega_{c,Z}}$")
#ax2.plot(theta, criteria_lcfs_good, lw=lw, color="tab:red")
ax2.plot(theta, criteria_length_lcfs, lw=lw+1, color="k")
ax2.plot(theta, criteria_length_lcfs, lw=lw, color="tab:brown", label=r"$\mathdefault{L_E / \rho_Z}$")
ax2.set_ylim([0.1, 10000])
ax2.set_yscale("log")
ax2.set_xlabel("Poloidal Angle", fontsize=fontsize)
ax2.set_ylabel("Guiding Center Validity", fontsize=fontsize)
ax2.legend(fontsize=fontsize-4, loc="lower right")

"""
ax3.scatter(char_t_lcfs * omega_c_lcfs, diff / v_drift_total)
ax3.set_xlabel("ratio")
ax3.set_ylabel("diff")

nX, binsX, patches = ax2.hist(criteria_X, bins=bins, density=True, 
	linewidth=1, edgecolor="k")
for patch, color in zip(patches, bar_colors):
	patch.set_facecolor(color)
ax2.set_xscale("log")
ax2.axvline(0.1, color="k", lw=3)
ax2.set_xlabel(r"$\mathdefault{\tau_{E_X} / \tau_{\omega_c}}$", 
	fontsize=fontsize)
ax2.text(0.75, 0.8, "{:.0f}%".format(frac_above(0.1, nX, binsX) * 100), 
	transform=ax2.transAxes, fontsize=fontsize, color="tab:red")

nY, binsY, patches = ax3.hist(criteria_Y, bins=bins, density=True, 
	linewidth=1, edgecolor="k")
for patch, color in zip(patches, bar_colors):
	patch.set_facecolor(color)
ax3.set_xscale("log")
ax3.axvline(0.1, color="k", lw=3)
ax3.set_xlabel(r"$\mathdefault{\tau_{E_Y} / \tau_{\omega_c}}$", 
	fontsize=fontsize)
ax3.text(0.75, 0.8, "{:.0f}%".format(frac_above(0.1, nY, binsY) * 100), 
	transform=ax3.transAxes, fontsize=fontsize, color="tab:red")

nZ, binsZ, patches = ax4.hist(criteria_Z, bins=bins, density=True, 
	linewidth=1, edgecolor="k")
for patch, color in zip(patches, bar_colors):
	patch.set_facecolor(color)
ax4.set_xscale("log")
ax4.axvline(0.1, color="k", lw=3)
ax4.set_xlabel(r"$\mathdefault{\tau_{E_Z} / \tau_{\omega_c}}$", 
	fontsize=fontsize)
ax4.text(0.75, 0.8, "{:.0f}%".format(frac_above(0.1, nZ, binsZ) * 100), 
	transform=ax4.transAxes, fontsize=fontsize, color="tab:red")
"""
fig.tight_layout()
fig.show()

# Plot of just Flan results of velocity
fig, ax1 = plt.subplots(figsize=(6, 5))

#for ax in [ax1, ax2, ax3, ax4]:
for ax in [ax1]:
	ax.tick_params(axis="both", which="both", labelsize=fontsize-2)

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

ax1.axhline(0.0, color="k")
ax1.plot(theta, v_perp_actual, lw=lw+3, color="k")
ax1.plot(theta, v_perp_actual, label="Flan", lw=lw+2, color="tab:cyan")
ax1.set_ylim([-10000, 10000])
#ax1.legend(fontsize=fontsize-6, ncols=2, loc="lower center")
ax1.set_xlabel("Poloidal Angle", fontsize=fontsize)
ax1.set_ylabel("Speed Leaving Core (m/s)", fontsize=fontsize) 

fig.tight_layout()
fig.show()
