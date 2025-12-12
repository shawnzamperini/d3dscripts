import get_lp
import numpy as np
import pickle
import matplotlib.pyplot as plt
import warnings
from scipy.signal import medfilt


# A list of key shots and the [tmin, tmax] to load data for
shots = {

	# LSN H-Mode
	200083: {"tminmax": [2000, 4000], 
		"exclude":[],
		"gfile_name": "200083_2400"},

	# Super-H
	203188: {"tminmax": [3000, 3300], 
		"exclude":[],
		"gfile_name": "203188_3000"},

	# LSN H-mode w/ baffle probes
	203720: {"tminmax": [2000, 4000], 
		"exclude":[],
		"gfile_name": "203720_3500"},

	# LSN L-mode w/ baffle probes
	203795: {"tminmax": [2000, 4000], 
		"exclude":[],
		"gfile_name": "203795_3500"},

	# QH-Mode
	201661: {"tminmax": [3000, 5000], 
		"exclude":[],
		"gfile_name": "201661_3000"},

	# NT
	193765: {"tminmax": [5200, 5400], 
		"exclude":[],
		"gfile_name": "193765_5500"},

	# High beta-p
	170361: {"tminmax": [4000, 5000], 
		"exclude":[],
		"gfile_name": "170361_4300"},

	# Helicon probes were on
	203665: {"tminmax": [2000, 4500], 
		"exclude":[],
		"gfile_name": "203665_3500"}
	}

# Time window for median filter
med_twin = 50  # ms

# Sputtering ion
sput_ion = "neon"
sput_ion_ch = 3

# Ti = fact * Te
ti_fact = 1.0

# Eckstein yield calculation for W
def Y_W(Ti, ion_ch, ion, debug=False):
	"""
	Fitted yield from Eckstein, assuming normal incidence. User only needs to
	enter the incident ion energy, E0, and what the ion is, either "deuterium"
	or "carbon".
		E0:
			Incident ion energy (eV).
		m1, m2:
			Mass of projectile and target atom, respectively.
		Z1, Z2:
			Atomic number of projectile and target atom, respectively.
		Eth, q, mu, lamb:
			Fitting parameters determined from Eckstein.
			Eth: Threshold energy for sputtering (eV).
			q: Related to absolute yield.
			lamb: Related to onset of decrease of yield at low energies towards
				  the threshold.
			mu: Describes strength of this decrease.
	"""

	# Make sure in lowercase.
	ion = ion.lower()

	# Impact energy is 5*Ti
	E0 = 5 * ion_ch * Ti

	# Masses in amu since all it uses it the ratio between them.
	m2 = 183.84
	Z2 = 74
	elec_sq = 1.44	# In eV*nm

	# Select the relevant fitting parameters.
	if ion == "deuterium":
		m1	 = 2.014
		Z1	 = 1
		q	 = 0.0183
		mu	 = 1.4410
		lamb = 0.3583
		Eth  = 228.84
	elif ion == "tungsten":
		m1	 = 183.84
		Z1	 = 74
		q	 = 18.6006
		mu	 = 3.1273
		lamb = 2.2697
		Eth  = 24.9885
	elif ion == "carbon":
		# Carbon parameters are not readily avaible, so they were determined
		# using the function below (i.e. fitting to TRIM data).
		m1	 = 12.01
		Z1	 = 6
		q	 = 5.4744
		mu	 = 2.0321
		lamb = 20.5545
		Eth  = 30.0
	elif ion == "neon":
		m1 = 20.18
		Z1 = 10
		q = 2.552
		mu = 1.9534
		lamb = 0.0828
		Eth = 38.6389
		#epsilon = 1.19E+05
		#Esb = 8.68
		
	else:
		print("Error: Ion must be one of:\n  Deuterium\n  Tungsten\n  Carbon")
		return None


	def lindhard():
		""" The Lindhard screening length (nm)."""
		return (9 * np.pi**2 / 128)**(1/3) * 0.0529177 * (Z1**(2/3) + Z2**(2/3))**(-1/2)

	def reduced_energy():
		""" The reduced energy (unitless). """
		return E0 * m1/(m1+m2) * lindhard() / (Z1 * Z2 * elec_sq)

	def nuclear_stopping():
		""" The nuclear stopping power with KrC (WHB) potential (). """
		return 0.5 * np.log(1 + 1.2288 * reduced_energy()) / (reduced_energy() +
			   0.1728 * np.sqrt(reduced_energy()) + 0.008 * reduced_energy()**0.1504)

	if debug:
		for index in range(0, len(E0)):
			print("Energy = {} eV".format(E0[index]))
			print("  Lindhard screening length = {:.4e}".format(lindhard()[index]))
			print("  Reduced energy =			 {:.4e}".format(reduced_energy()[index]))
			print("  Nuclear stopping power =	 {:.4e}".format(nuclear_stopping()[index]))
			print("")

	# Return the yield using the above functions.
	# Weird warning only involving np.floats, but doesn't seem to cause problems so ignore it.
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		ans = q * nuclear_stopping() * (E0 / Eth - 1)**mu / (lamb + (E0 / Eth - 1)**mu)
		if isinstance(ans, complex):
			return 0
		else:
			return ans

def mask_data(data_dict, mask):
	output = {}
	for k, v in data_dict.items():
		output[k] = np.array(v)[mask]
	return output

# Figure to plot a particular probe for inspection
probe_plot = "probe 37"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))

shot_data = {}
for shot, aux in shots.items():
	print("-- {} --".format(shot))
	tminmax = aux["tminmax"]
	exclude = aux["exclude"]
	
	# Load gfile for magnetic data
	with open("/home/zamp/data/{}.pickle".format(aux["gfile_name"]), "rb") as f:
		gfile = pickle.load(f)

	# Pull out arrays for conveience
	Bp = gfile["Bp"]
	Bz = gfile["Bz"]
	Br = gfile["Br"]
	Bt = gfile["Bt"]
	psin = gfile["PSIRZ_NORM"]
	B = np.sqrt(np.square(Br) + np.square(Bz) + np.square(Bt))

	# Load LP data
	lps = get_lp.get_dict_of_lps(shot, tunnel=False)

	# Go through each probe
	med_te = []
	med_ne = []
	med_jsat = []
	avg_t = []
	b = []
	angs = []
	r = []
	z = []
	wall_bp = []
	wall_b = []
	wall_psin = []
	y = []
	labels = []
	for pname, lp in lps.items():
		
		# If probe excluded, skip
		if lp["label"].strip() in exclude: continue

		# Pull out data we care about
		t = lp["time"]
		te = lp["temp"]
		ne = lp["dens"]
		ang = lp["angle"]
		isat = lp["isat"]
		jsat = lp["jsat"]
		psin = lp["psin"]
		csq = lp["csq"]
		label = np.full(len(t), lp["label"].strip())
		data_dict = {"t":t, "te":te, "ne":ne, "ang":ang, "isat":isat, 
			"jsat":jsat, "psin":psin, "csq":csq, "label":label}

		if pname == probe_plot:
			ax1.plot(data_dict["t"], data_dict["te"])
			ax2.plot(data_dict["t"], data_dict["jsat"])

		# Limit to tmin, tmax
		mask = np.logical_and(data_dict["t"] > tminmax[0], 
			data_dict["t"] < tminmax[1])
		data_dict = mask_data(data_dict, mask)

		# Only use values where jsat, Te and psin are not equal to zero
		mask = data_dict["jsat"] != 0
		data_dict = mask_data(data_dict, mask)
		mask = data_dict["te"] != 0
		data_dict = mask_data(data_dict, mask)
		mask = data_dict["psin"] != 0
		data_dict = mask_data(data_dict, mask)

		# Values where isat is not -0.05 (an error code apparently)
		mask = data_dict["isat"] != -0.05
		data_dict = mask_data(data_dict, mask)

		# Values where Te is not equal to 0.1
		mask = data_dict["te"] != 0.01
		data_dict = mask_data(data_dict, mask)

		# Values where Te is not greater than 200
		mask = data_dict["te"] < 200
		data_dict = mask_data(data_dict, mask)

		# Values where jsat > 0.001 A/cm2
		mask = data_dict["jsat"] > 0.001
		data_dict = mask_data(data_dict, mask)

		# Values where isat <= 4A
		mask = data_dict["isat"] <= 4.0
		data_dict = mask_data(data_dict, mask)

		# Values where chisquared is not 99999
		mask = data_dict["csq"] != 99999
		data_dict = mask_data(data_dict, mask)

		if pname == probe_plot:
			ax1.plot(data_dict["t"], data_dict["te"])
			ax2.plot(data_dict["t"], data_dict["jsat"])

		# Now apply a median filter to the data
		medfilt_window = 51
		data_dict["te"] = medfilt(data_dict["te"], medfilt_window)
		data_dict["ne"] = medfilt(data_dict["ne"], medfilt_window)
		data_dict["jsat"] = medfilt(data_dict["jsat"], medfilt_window)
		
		if pname == probe_plot:
			ax1.plot(data_dict["t"], data_dict["te"])
			ax2.plot(data_dict["t"], data_dict["ne"])

		# Use gfile to find the B strength
		# Calculate distance to each gfile data point and get index of minimum
		dist = np.sqrt(np.square(gfile["R"] - lp["rprobe"]) 
			+ np.square(gfile["Z"] - lp["zprobe"]))
		min_index = np.unravel_index(np.argmin(dist), dist.shape)

		# Apply a median filter
		med_t = 0
		while med_t < t.max():

			# Mask within median window size
			mask = np.logical_and(data_dict["t"] >= med_t,
				data_dict["t"] < med_t + med_twin)
			if mask.sum() == 0: 
				med_t += med_twin
				continue
		
			# Calculate median plasma data
			#tmp_te = np.median(data_dict["te"][mask])
			#tmp_ne = np.median(data_dict["ne"][mask])
			tmp_te = np.mean(data_dict["te"][mask])
			tmp_ne = np.mean(data_dict["ne"][mask])
			tmp_jsat = np.mean(data_dict["jsat"][mask])
			med_te.append(tmp_te)
			med_ne.append(tmp_ne)
			med_jsat.append(tmp_jsat)

			if pname == probe_plot:
				tmp_t = data_dict["t"][mask].mean()
				ax1.scatter(tmp_t, tmp_te, color="k", zorder=15)
				ax2.scatter(tmp_t, tmp_jsat, color="k", zorder=15)

			# Calculate sputtering yield
			tmp_y = Y_W(tmp_te * ti_fact, sput_ion_ch, sput_ion)
			y.append(tmp_y)
			#print("{:6.2f} {:8.2e}".format(tmp_te, tmp_y))

			# Logical to append average t to each data point
			avg_t.append(data_dict["t"][mask].mean())

			# Add R, Z, average angle data
			r.append(lp["rprobe"])
			z.append(lp["zprobe"])
			angs.append(data_dict["ang"][mask].mean())

			# Get Bp and B
			wall_bp.append(Bp[min_index])
			wall_b.append(B[min_index])

			# Get psin
			wall_psin.append(psin[min_index])

			# Append label so every data point can be identified.
			labels.append(data_dict["label"][0])

			# Increment for next time window
			med_t += med_twin

	# Store in shot_data for post-processing
	shot_data[shot] = {"t":avg_t, "te":med_te, "ne":med_ne, "r":r, "z":z, 
		"angs":angs, "wall_bp":wall_bp, "wall_b":wall_b, "wall_psin":wall_psin,
		"y":y, "jsat":med_jsat, "labels":labels}

fig.tight_layout()
fig.show()

# Additional arrays with the data consolidated
all_te = []
all_ne = []
all_jsat = []
all_r = []
all_z = []
all_y = []
all_labels = []
for shot, data in shot_data.items():
	all_te = np.append(all_te, data["te"])
	all_ne = np.append(all_ne, data["ne"])
	all_jsat = np.append(all_jsat, data["jsat"])
	all_r = np.append(all_r, data["r"])
	all_z = np.append(all_z, data["z"])
	all_y = np.append(all_y, data["y"])
	all_labels = np.append(all_labels, data["labels"])
all_te = np.array(all_te).flatten()
all_ne = np.array(all_ne).flatten()
all_jsat = np.array(all_jsat).flatten()
all_r = np.array(all_r).flatten()
all_z = np.array(all_z).flatten()
all_y = np.array(all_y).flatten()
all_labels = np.array(all_labels).flatten()

# nans can be zeros
all_y[np.isnan(all_y)] = 0.0

# Calculate erosion flux = Y * jsat * factor
all_ero = all_y * all_jsat * 6.2415e18	# eroded particles/s/cm2

# Ignore all data below a threshold
mask = all_ero > np.max(all_ero) * 0.0001
all_te = all_te[mask]
all_ne = all_ne[mask]
all_jsat = all_jsat[mask]
all_r = all_r[mask]
all_z = all_z[mask]
all_y = all_y[mask]
all_ero = all_ero[mask]
all_labels = all_labels[mask]

# Sort each, creating a matching R, Z array for the plot
te_idx = np.argsort(all_ero)[::-1]
all_te = all_te[te_idx]
all_te_r = all_r[te_idx]
all_te_z = all_z[te_idx]
all_te_y = all_y[te_idx]
all_te_ero = all_ero[te_idx]
all_te_labels = all_labels[te_idx]

# Load a gfile to grab the wall coordinates
with open("/home/zamp/data/203188_3000.pickle", "rb") as f:
	gfile = pickle.load(f)

# Figure with wall coordinates
fig, ax1 = plt.subplots(figsize=(7, 10))
ax1.set_aspect("equal")
ax1.plot(gfile["RLIM"], gfile["ZLIM"], color="k", lw=3, zorder=5)

# Normalize size between smin, smax
def normalize(data, new_min, new_max, old_min=None, old_max=None):
	if old_min is None:
		old_min = np.min(data)
	if old_max is None:
		old_max = np.max(data)
	return new_min + ((data - old_min) * (new_max - new_min) / (old_max - old_min))

def log_normalize(data, new_min=0, new_max=1, old_min_mod=1.0):
	"""Normalize an array of values using logarithmic scaling between a set range."""
	if old_min_mod > 1.0:
		print("Warning! old_min_mod should be <0, otherwise you will lose data points.")
	data = np.array(data)
	data = np.where(data > 0, np.log(data), np.nan)  # Apply log transformation to positive values
	log_min, log_max = np.nanmin(data * old_min_mod), np.nanmax(data)  # Get min/max values ignoring NaNs
	normalized_data = new_min + ((data - log_min) * (new_max - new_min) / (log_max - log_min))	# Scale to range
	return normalized_data


smin = 5
smax = 750
#s = normalize(all_te_ero, smin, smax)
s = log_normalize(all_te_ero, smin, smax)
#s = log_normalize(all_ero, smin, smax)
#s[s<0] = smin

# Plot circles normalized to the maximum data value
ax1.scatter(all_te_r, all_te_z, s=s, color="tab:red", edgecolors="k", zorder=15)
#ax1.scatter(all_r, all_z, s=s, color="tab:red", edgecolors="k", zorder=15)

fig.tight_layout()
fig.show()


# Another way to show this data, use a colormap to colorcode each wall segment
# based on the nearest data point. 
rlim = gfile["RLIM"]
zlim = gfile["ZLIM"]

# Break each segment up into more discrete segments. This is primarily for the
# large straight regions like the inner wall
nbreaks = 10
rlim_hires = []
zlim_hires = []
for i in range(len(rlim) - 1):
	
	# Point-slope formula. Account for vertical segments.
	if (np.abs(rlim[i+1] - rlim[i]) < 1e-10):
		for j in range(nbreaks):
			rlim_hires.append(rlim[i])
			zlim_hires.append(zlim[i])
		continue
	else:
		m = (zlim[i+1] - zlim[i]) / (rlim[i+1] - rlim[i])

	# Incrementally add newer points
	for j in range(nbreaks):
		newr = rlim[i] + (rlim[i+1] - rlim[i]) / nbreaks * j
		newz = zlim[i] + m * (newr - rlim[i])
		rlim_hires.append(newr)
		zlim_hires.append(newz)

# Don't forget to append the last coordinate.
rlim_hires.append(rlim[-1])
zlim_hires.append(zlim[-1])

# Calculate the midpoint of each segment
rmid_hires = np.array([(rlim_hires[i] + rlim_hires[i+1]) / 2.0 for i in range(len(rlim_hires)-1)])
zmid_hires = np.array([(zlim_hires[i] + zlim_hires[i+1]) / 2.0 for i in range(len(zlim_hires)-1)])

# Only consider maximum values for each LP location, identified by label
max_eros = {}
for label in np.unique(all_te_labels):
	mask = all_te_labels == label
	max_eros[label] = {"r":all_te_r[mask][0], "z":all_te_z[mask][0], 
		"ero":all_te_ero[mask].max()}

# Find nearest LP data point
hires_ero = []
for i in range(len(rmid_hires)):

	# Find nearest probe to the midpoint value
	nearest_label = ""
	mindist = 1e12
	for k, d in max_eros.items():
		dist = np.sqrt(np.square(rmid_hires[i] - d["r"]) 
			+ np.square(zmid_hires[i] - d["z"]))	
		if dist < mindist:
			nearest_label = k
			mindist = dist

	# Append erosion value	
	hires_ero.append(max_eros[nearest_label]["ero"])

# Now create a plot where the color of the line segment corresponds to the
# erosion value
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LogNorm

# Create list of segments
segments = np.array([[[rlim_hires[i], zlim_hires[i]], 
	[rlim_hires[i+1], zlim_hires[i+1]]] for i in range(len(rlim_hires)-1)])

# Colormapping based on erosion value
#norm = Normalize(vmin = min(hires_ero), vmax = max(hires_ero))
norm1 = LogNorm(vmin = min(hires_ero), vmax = max(hires_ero))
#norm2 = LogNorm(vmin = min(hires_ero) / 10, vmax = max(hires_ero) / 10)
cmap = plt.cm.Reds_r
colors1 = cmap(norm1(hires_ero[:-1]))
#colors2 = cmap(norm2(hires_ero[:-1]))

# Create a LineCollection
lc1 = LineCollection(segments, colors=colors1, linewidths=6)
#lc2 = LineCollection(segments, colors=colors2, linewidths=6)

# Make plot
fig, ax1 = plt.subplots(figsize=(7, 10))
ax1.add_collection(lc1)
#ax2.add_collection(lc2)
for ax in [ax1]:
	ax.set_facecolor("k")
	ax.set_xlim([0.95, 2.5])
	ax.set_ylim([-1.5, 1.5])
	ax.set_aspect("equal")

# Add colorbar
sm1 = plt.cm.ScalarMappable(norm=norm1, cmap=cmap)
#sm2 = plt.cm.ScalarMappable(norm=norm2, cmap=cmap)
sm1.set_array([])
#sm2.set_array([])
cbar1 = fig.colorbar(sm1, ax=ax1)
#cbar2 = fig.colorbar(sm2, ax=ax2)
cbar1.ax.set_ylabel("W Erosion (W/m2/s)", fontsize=16)
#cbar2.ax.set_ylabel("W Erosion (W/m2/s)", fontsize=16)
cbar1.ax.tick_params(labelsize=14)
#cbar2.ax.tick_params(labelsize=14)

plt.show()
