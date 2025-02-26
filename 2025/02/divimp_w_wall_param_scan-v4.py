import oedge_plots
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression


# Figure of merit. This is what is used to compare the different scenarios
# to each other. Options are: avg_conc, avg_nw
fom = "avg_conc"

# Choose the ring just inside the separatrix for analysis
ring = 15

# Create figure for the plotting below
fontsize = 16
fig, ax = plt.subplots()

# The parameters that were scanned
dperps = [0.593, 3.900, 2.760, 4.864, 4.446, 2.840, 1.248, 1.110, 2.257, 3.557, 3.647, 0.304, 1.463, 2.778, 0.632, 4.843, 2.137, 0.420, 3.211, 3.632, 3.542, 2.489, 1.491, 0.546, 3.899, 1.548, 2.292, 1.631, 3.412, 3.136, 0.180, 0.542, 0.225, 3.196, 1.643, 2.967, 2.747, 4.074, 3.153, 0.805, 0.161, 3.167, 3.359, 4.068, 0.596, 4.870, 0.321, 0.313, 3.232, 4.666, 2.524, 0.424, 2.559, 3.792, 0.333, 1.986, 1.168, 3.398, 3.989, 4.612, 1.898, 0.337, 2.796, 3.017, 3.093, 2.152, 1.658, 0.263, 3.612, 4.772, 2.978, 0.577, 1.129, 2.490, 4.321, 0.418, 2.435, 0.294, 2.519, 1.106, 0.468, 1.033, 1.388, 2.710, 3.658, 2.158, 4.680, 3.592, 0.981, 2.466, 4.886, 3.915, 3.525, 2.242, 1.293, 2.358, 2.508, 0.929, 4.337, 2.856, 4.326, 3.750, 3.485, 2.005, 0.202, 2.928, 2.000, 1.841, 3.044, 2.077, 0.152, 4.461, 2.313, 3.070, 1.236, 4.346, 4.429, 2.664, 3.647, 1.896, 3.346, 0.833, 1.122, 4.822, 0.935, 1.652, 2.874, 1.321, 2.476, 3.171, 4.942, 4.825, 1.125, 4.089, 0.444, 2.431, 4.286, 2.526, 2.869, 3.694, 2.741, 4.445, 3.959, 1.231, 0.982, 3.133, 4.732, 0.249, 3.431, 0.712, 2.243, 4.973, 0.371, 2.579, 0.227, 2.143, 2.638, 4.241, 4.879, 4.479, 4.009, 1.638, 2.543, 3.332, 2.379, 1.678, 0.491, 4.967, 3.058, 1.873, 4.317, 3.946, 0.676, 3.191, 2.669, 0.242, 4.122, 4.279, 4.161, 4.527, 0.135, 1.068, 2.904, 0.330, 3.608, 4.906, 3.841, 1.935, 0.528, 4.130, 2.476, 4.951, 2.709, 4.759, 1.093, 2.229, 2.631, 3.023, 4.911, 2.986, 0.702, 4.361, 1.296, 0.400, 4.500, 3.878, 3.916, 4.304, 4.089, 1.163, 4.276, 1.518, 1.380, 3.499, 4.695, 0.934, 3.301, 1.060, 1.406, 0.453, 1.044, 4.757, 3.388, 3.802, 3.041, 4.460, 4.918, 3.379, 0.618, 1.796, 2.466, 3.443, 4.917, 0.274, 2.126, 3.469, 2.659, 1.244, 3.310, 0.601, 0.191, 3.058, 2.275, 1.058, 1.551, 2.968, 3.077, 2.556, 2.092, 1.256]
ne_charges = [5, 2, 6, 6, 9, 3, 3, 9, 6, 5, 2, 9, 6, 6, 6, 3, 4, 2, 5, 8, 7, 8, 4, 6, 2, 2, 3, 7, 7, 4, 4, 3, 9, 4, 8, 9, 7, 5, 2, 7, 3, 5, 9, 8, 4, 7, 4, 6, 7, 9, 6, 2, 3, 4, 9, 8, 9, 6, 4, 5, 8, 9, 5, 6, 7, 4, 6, 6, 3, 5, 4, 9, 6, 8, 8, 9, 2, 5, 6, 2, 7, 7, 9, 2, 7, 8, 3, 9, 9, 2, 2, 3, 3, 8, 2, 7, 8, 5, 4, 6, 7, 3, 7, 2, 9, 4, 5, 8, 4, 6, 6, 5, 8, 7, 9, 8, 3, 7, 6, 4, 3, 4, 3, 4, 3, 6, 8, 5, 7, 6, 9, 6, 4, 6, 2, 5, 6, 2, 2, 3, 2, 7, 4, 4, 8, 5, 5, 5, 8, 2, 3, 4, 6, 2, 3, 6, 3, 8, 6, 9, 4, 6, 9, 9, 5, 6, 8, 4, 9, 8, 8, 8, 8, 6, 4, 9, 9, 2, 5, 7, 7, 2, 9, 7, 7, 3, 6, 3, 3, 4, 7, 7, 4, 6, 6, 9, 7, 2, 5, 9, 7, 3, 9, 2, 4, 6, 9, 5, 5, 5, 6, 3, 5, 7, 9, 5, 9, 9, 4, 6, 2, 6, 9, 2, 6, 8, 9, 8, 4, 2, 2, 5, 5, 7, 8, 3, 8, 2, 5, 5, 9, 5, 3, 9, 9, 3, 9, 4, 8, 5]
ne_fluxes = [0.006, 0.058, 0.040, 0.093, 0.059, 0.017, 0.039, 0.077, 0.097, 0.092, 0.073, 0.071, 0.008, 0.015, 0.014, 0.072, 0.054, 0.095, 0.009, 0.048, 0.065, 0.039, 0.062, 0.006, 0.056, 0.009, 0.026, 0.060, 0.049, 0.066, 0.016, 0.097, 0.087, 0.094, 0.058, 0.052, 0.099, 0.083, 0.047, 0.073, 0.062, 0.047, 0.046, 0.074, 0.022, 0.025, 0.041, 0.075, 0.054, 0.056, 0.013, 0.094, 0.061, 0.067, 0.031, 0.078, 0.096, 0.086, 0.049, 0.051, 0.046, 0.052, 0.037, 0.034, 0.088, 0.081, 0.041, 0.009, 0.026, 0.014, 0.041, 0.041, 0.015, 0.089, 0.089, 0.014, 0.082, 0.067, 0.096, 0.089, 0.055, 0.094, 0.066, 0.025, 0.019, 0.099, 0.029, 0.031, 0.092, 0.040, 0.038, 0.024, 0.043, 0.012, 0.057, 0.047, 0.097, 0.071, 0.034, 0.088, 0.049, 0.025, 0.098, 0.033, 0.078, 0.057, 0.010, 0.091, 0.084, 0.093, 0.056, 0.099, 0.040, 0.064, 0.064, 0.019, 0.024, 0.039, 0.064, 0.011, 0.099, 0.085, 0.059, 0.099, 0.023, 0.086, 0.014, 0.084, 0.037, 0.040, 0.061, 0.038, 0.071, 0.022, 0.093, 0.089, 0.048, 0.035, 0.088, 0.052, 0.080, 0.080, 0.091, 0.023, 0.044, 0.010, 0.039, 0.051, 0.047, 0.055, 0.011, 0.032, 0.098, 0.079, 0.020, 0.032, 0.032, 0.072, 0.008, 0.043, 0.055, 0.088, 0.056, 0.068, 0.073, 0.088, 0.099, 0.041, 0.013, 0.039, 0.013, 0.017, 0.022, 0.039, 0.044, 0.077, 0.089, 0.059, 0.007, 0.071, 0.074, 0.025, 0.025, 0.089, 0.063, 0.098, 0.076, 0.014, 0.054, 0.062, 0.026, 0.084, 0.017, 0.063, 0.049, 0.082, 0.033, 0.091, 0.022, 0.009, 0.034, 0.026, 0.039, 0.059, 0.060, 0.083, 0.012, 0.096, 0.094, 0.085, 0.076, 0.067, 0.092, 0.085, 0.069, 0.090, 0.069, 0.010, 0.100, 0.086, 0.075, 0.020, 0.027, 0.059, 0.064, 0.028, 0.097, 0.007, 0.052, 0.070, 0.093, 0.009, 0.057, 0.099, 0.072, 0.031, 0.019, 0.075, 0.086, 0.091, 0.083, 0.075, 0.058, 0.020, 0.096, 0.083, 0.013, 0.076, 0.086, 0.026]
wall_tes = [7.919, 2.710, 6.509, 3.165, 8.277, 2.357, 4.670, 1.124, 4.547, 3.937, 4.776, 9.293, 3.749, 2.715, 6.454, 8.780, 2.530, 1.382, 5.444, 1.014, 5.726, 1.319, 6.194, 8.345, 3.065, 9.447, 5.287, 7.791, 1.920, 4.986, 9.833, 3.310, 8.923, 7.718, 4.900, 9.577, 3.000, 9.363, 4.020, 7.831, 2.914, 8.602, 4.566, 8.326, 4.962, 6.475, 1.269, 3.435, 8.700, 3.253, 5.487, 5.374, 5.497, 3.494, 1.829, 4.647, 4.061, 8.162, 1.233, 9.313, 9.395, 6.377, 4.164, 3.004, 3.805, 8.277, 5.682, 8.690, 4.158, 1.193, 2.295, 1.953, 7.881, 8.623, 8.491, 3.541, 4.108, 8.253, 3.526, 5.314, 8.944, 8.816, 9.418, 7.300, 4.798, 9.485, 9.394, 9.017, 1.741, 5.054, 1.583, 6.639, 8.949, 6.730, 6.062, 4.992, 8.371, 5.368, 1.252, 5.919, 6.271, 3.390, 4.601, 3.284, 4.967, 5.949, 3.694, 2.782, 2.884, 2.645, 7.404, 4.513, 5.444, 9.531, 3.556, 8.622, 1.798, 5.088, 3.601, 5.225, 9.574, 3.885, 2.111, 2.083, 6.475, 7.956, 6.616, 1.431, 9.636, 3.906, 8.215, 6.215, 8.673, 7.924, 3.868, 1.627, 4.271, 4.750, 3.935, 8.071, 4.194, 6.996, 5.739, 4.768, 7.721, 5.834, 5.083, 9.789, 6.080, 8.336, 3.304, 5.189, 4.335, 2.295, 2.685, 4.786, 2.719, 2.501, 7.969, 9.166, 2.891, 3.798, 5.438, 7.591, 4.140, 5.581, 9.256, 6.712, 3.048, 2.309, 3.992, 8.678, 5.231, 8.914, 3.900, 1.626, 8.152, 6.318, 8.514, 5.088, 5.622, 7.175, 6.337, 5.162, 6.501, 8.932, 4.615, 6.797, 7.233, 8.375, 2.044, 9.597, 5.413, 6.962, 4.374, 9.059, 5.507, 6.965, 4.931, 3.240, 2.161, 2.712, 4.297, 8.410, 8.269, 6.760, 4.568, 5.551, 9.706, 1.678, 1.412, 1.969, 2.962, 5.281, 5.672, 7.488, 7.424, 7.733, 2.470, 2.157, 7.871, 7.611, 3.692, 8.492, 9.759, 7.008, 1.613, 9.956, 3.722, 6.613, 9.616, 1.806, 3.441, 5.053, 4.489, 5.161, 4.391, 5.684, 8.652, 7.638, 7.522, 1.128, 5.893, 5.034, 2.218, 8.625, 2.871, 5.392, 8.967, 4.917]
machs = [-0.087, -0.455, -0.117, -0.396, -0.117, -0.157, -0.040, -0.251, -0.483, -0.370, -0.099, -0.238, -0.296, -0.464, -0.413, -0.194, -0.348, -0.309, -0.128, -0.177, -0.113, -0.187, -0.304, -0.207, -0.366, -0.480, -0.016, -0.333, -0.448, -0.101, -0.352, -0.258, -0.470, -0.287, -0.196, -0.050, -0.458, -0.124, -0.033, -0.470, -0.358, -0.240, -0.290, -0.193, -0.455, -0.206, -0.256, -0.462, -0.390, -0.284, -0.330, -0.438, -0.019, -0.152, -0.026, -0.358, -0.068, -0.152, -0.124, -0.335, -0.301, -0.440, -0.343, -0.022, -0.308, -0.077, -0.374, -0.016, -0.033, -0.498, -0.076, -0.470, -0.232, -0.422, -0.298, -0.095, -0.037, -0.004, -0.226, -0.083, -0.362, -0.002, -0.462, -0.238, -0.082, -0.276, -0.328, -0.444, -0.002, -0.192, -0.138, -0.329, -0.186, -0.407, -0.242, -0.287, -0.405, -0.012, -0.285, -0.052, -0.305, -0.044, -0.313, -0.335, -0.312, -0.282, -0.058, -0.071, -0.483, -0.056, -0.365, -0.437, -0.266, -0.432, -0.248, -0.088, -0.277, -0.081, -0.191, -0.254, -0.408, -0.248, -0.007, -0.345, -0.060, -0.330, -0.164, -0.240, -0.361, -0.161, -0.461, -0.243, -0.081, -0.245, -0.254, -0.169, -0.022, -0.152, -0.356, -0.164, -0.234, -0.312, -0.298, -0.104, -0.287, -0.051, -0.354, -0.228, -0.117, -0.471, -0.042, -0.130, -0.032, -0.384, -0.081, -0.374, -0.106, -0.446, -0.137, -0.400, -0.016, -0.185, -0.382, -0.395, -0.100, -0.318, -0.174, -0.014, -0.252, -0.388, -0.073, -0.381, -0.251, -0.405, -0.097, -0.423, -0.069, -0.313, -0.145, -0.228, -0.114, -0.128, -0.131, -0.258, -0.383, -0.411, -0.106, -0.243, -0.256, -0.351, -0.353, -0.473, -0.382, -0.288, -0.253, -0.322, -0.344, -0.229, -0.437, -0.104, -0.047, -0.382, -0.122, -0.236, -0.413, -0.384, -0.264, -0.170, -0.111, -0.067, -0.046, -0.404, -0.191, -0.162, -0.050, -0.221, -0.145, -0.028, -0.143, -0.486, -0.219, -0.046, -0.379, -0.474, -0.315, -0.279, -0.031, -0.192, -0.482, -0.416, -0.414, -0.289, -0.085, -0.017, -0.398, -0.481, -0.276, -0.017, -0.108, -0.081, -0.339, -0.363, -0.145, -0.254, -0.361, -0.467, -0.055, -0.444, -0.143, -0.412]

# Zip together so we have a record of each case's setting
xvalues = np.array(list(zip(dperps, ne_charges, ne_fluxes, wall_tes, machs)))
#xvalues = xvalues[:50]

# Make Mach values positive, arbitrary
xvalues[:,4] *= -1

# Normalize the xvalues between 0-1
#for i in range(0, xvalues.shape[1]):
#	xvalues[:,i] = (xvalues[:,i] - xvalues[:,i].min()) / \
#		(xvalues[:,i].max() - xvalues[:,i].min())

# Normalizing data
scaler = StandardScaler()
xvalues_norm = scaler.fit_transform(xvalues)

# Some quick 2D plots that are used in the presentation
ncpath = "/home/zamp/oedge_files/d3d-w-wall-param-scan-v4-all-1.nc"
op = oedge_plots.OedgePlots(ncpath)
op.plot_contour_polygon("KTEBS", normtype="log", cbar_label="Te (eV)", 
	cmap="inferno", vmin=1)
op.plot_contour_polygon("DDLIMS", normtype="log", cbar_label="W Density (m-3)",
	charge="all", cmap="inferno", scaling=op.absfac, vmin=1e16, vmax=1e19)

fom_vals_dict = {}
fom_stds_dict = {}
#sput_regions = ["all", "ub", "ow", "lb", "sh", "iw"]
sput_regions = ["all"]
for sput_region in sput_regions:

	# Lists to hold our figures of merit as we extract them from each case
	fom_vals = []
	fom_stds = []

	for i in range(0, xvalues.shape[0]):

		# Load case
		ncpath = ("/home/zamp/oedge_files/d3d-w-wall-param-scan-v4-{}-{}.nc"
			.format(sput_region, i+1))
		print(ncpath)
		op = oedge_plots.OedgePlots(ncpath)

		# Load needed density to derive figure of merit at desired ring
		s, nw_s = op.along_ring(ring, "DDLIMS", charge="all", 
			scaling=op.absfac, plot_it=False)
		s, ne_s = op.along_ring(ring, "KNBS", plot_it=False)
		cw_s = nw_s / (nw_s + ne_s)

		# Figure of merit: Average W concentration
		if fom == "avg_conc":
			ax.plot(s, cw_s, lw=4, alpha=0.3, color="k")
			ax.plot(s, cw_s, lw=3, alpha=0.3, color="tab:red")

			fom_vals.append(cw_s.mean())
			fom_stds.append(cw_s.std())
		
		# Figure of merit: Average W density
		elif fom == "avg_nw":
			ax.plot(s, nw_s, lw=3, alpha=0.3, color="tab:red")

			fom_vals.append(nw_s.mean())
			fom_stds.append(nw_s.std())

	fom_vals_dict[sput_region] = np.array(fom_vals)
	fom_stds_dict[sput_region] = np.array(fom_stds)

#ax.set_yscale("log")
ax.set_xlabel("Distance from inner target (m)", fontsize=fontsize)

# Format y axis based off figure of merit.
if fom == "avg_conc":
	ax.set_ylabel("W Concentration", fontsize=fontsize)
	#ax.set_ylim([0, 1e-5])
elif fom == "avg_nw":
	ax.set_ylabel("W Density (m-3)", fontsize=fontsize)

#ax.set_yscale("log")
#ax.set_ylim([5e-6, 3e-2])
fig.tight_layout()
fig.show()


# Now that we have a large enough set of values, let's perform a multivariable
# curve fit to derive a scaling law
def model(x, a, b, c, d, e, f):
	
	# Extract each independent variable
	x1, x2, x3, x4, x5 = x

	# Simple power law for each variable
	return x1**a * x2**b * x3**c * x4**d * x5**e * f

fit_dict = {}
for sput_region in sput_regions:

	# Fit to our model function, save in dictionary
	popt, pcov = curve_fit(model, xvalues.T, fom_vals_dict[sput_region],
		bounds=((-10, -10, -10, -10, -10, -np.inf), 
		(10, 10, 10, 10, 10, np.inf)))
	fit_dict[sput_region] = model(xvalues.T, *popt)
	fit_dict[sput_region + "_popt"] = popt
	
	print("Fit Parameters ({})".format(sput_region))
	print("  dperp:     {:.2f}".format(popt[0]))
	print("  ne_charge: {:.2f}".format(popt[1]))
	print("  ne_flux:   {:.2f}".format(popt[2]))
	print("  wall_te:   {:.2f}".format(popt[3]))
	print("  mach:      {:.2f}".format(popt[4]))

# Can also fit a linear regression to the data
reg_dict = {}
for sput_region in sput_regions:
	linreg = LinearRegression().fit(xvalues_norm, fom_vals_dict[sput_region])
	reg_dict[sput_region] = linreg.predict(xvalues_norm)

# Used in making the x-axis label
if fom == "avg_conc":
	lhs = r"$\mathdefault{c_W}$"
elif fom == "avg_nw":
	lhs = r"$\mathdefault{n_W}$"

# Still growing this plot, but it's the figure of merit plot
fontsize = 18

# Make a plot for the two types of fits we are doing
for d in [fit_dict, reg_dict]:
	fig, axs = plt.subplots(2, 3, figsize=(14, 9))
	axs = axs.flatten()
	for i in range(len(sput_regions)):
		
		# Determines the limits of the axes
		lims = [1e-6, 1]

		# Plot data for this sputtering region
		key = sput_regions[i]
		axs[i].scatter(d[key], fom_vals_dict[key], s=75, 
			color="tab:red", edgecolors="k")
		axs[i].plot(lims, lims, lw=3, color="k", linestyle="--")

		# Create a text object to place on the plot of the fit
		if d == fit_dict:
			popt = fit_dict[key + "_popt"]
			text = lhs + r"$\mathdefault{\sim D_r^{" + "{:.2f}".format(popt[0]) \
				+ r"}\ Z_{Ne}^{" + "{:.2f}".format(popt[1]) \
				+ r"}\ f_{Ne}^{" + "{:.2f}".format(popt[2]) \
				+ r"}\ T_e^{" + "{:.2f}".format(popt[3]) \
				+ r"}\ M_{||}^{" + "{:.2f}".format(popt[4]) \
				+ r"}}$"
			axs[i].set_xlabel(text, fontsize=fontsize)
		elif d == reg_dict:

			# Extract coefficients and intercept
			c = linreg.coef_
			inter = linreg.intercept_

			# Build the linear regression equation string
			#terms = [f"{coeff:.3f} * x{i}" for i, coeff in enumerate(coefficients)]
			#equation = " + ".join(terms) + f" + {intercept:.3f}"

			text = lhs + " = {:.3f}".format(c[0]) + r"$\mathdefault{D_r" \
				+ " {:+.3f}".format(c[1]) + r"Z_{Ne}" \
				+ " {:+.3f}".format(c[2]) + r"f_{Ne}" \
				+ r"}$" + "\n" + r"$\mathdefault{" \
				+ " {:+.3f}".format(c[3]) + r"T_e" \
				+ " {:+.3f}".format(c[4]) + r"M_{||}" \
				+ " {:+.3f}".format(inter) + r"}$"
			axs[i].set_xlabel(text, fontsize=fontsize-2)

		# Plot formatting
		axs[i].set_title(key, fontsize=fontsize)
		axs[i].set_ylabel(lhs, fontsize=fontsize)
		axs[i].tick_params(axis="both", labelsize=fontsize-4)
		axs[i].grid(alpha=0.3)
		axs[i].set_xscale("log")
		axs[i].set_yscale("log")
		axs[i].set_xlim(lims)
		axs[i].set_ylim(lims)

	fig.tight_layout()
	fig.show()

# Need to run these simulations still
"""
# Now the task is to create a pie chart of the average contribution to the
# "all" figure of merit from each individual region.
region_fracs = []
region_fracs_dict = {"other": np.ones(len(dperps))}
for sput_region in sput_regions:
	if sput_region == "all": continue

	# Yay numpy arrays
	region_fracs_dict[sput_region] = fom_vals_dict[sput_region] / \
		fom_vals_dict["all"]
	region_fracs.append(region_fracs_dict[sput_region].mean())

	# Subtract from the other, whatever is left after the loop must come from
	# the regions that didn't get their own simulations
	region_fracs_dict["other"] -= region_fracs_dict[sput_region]

# Add on the other contribution, the convert to numpy array.
region_fracs.append(region_fracs_dict["other"].mean())
region_fracs = np.array(region_fracs)

pie_labels = ["Upper Baffle", "Outer Wall", "Lower Baffle", "Shelf", 
	"Inner Wall"]
fig, ax = plt.subplots(figsize=(10, 6))

# We don't include "other" because it is actually negative (about -0.005),
# which I interpret as a negligible contribution.
ax.pie(region_fracs[:-1], labels=pie_labels, autopct='%1.0f%%', 
	textprops={"fontsize": fontsize}, pctdistance=0.7, 
	wedgeprops={"edgecolor": "k", "linewidth": 2}, colors=["lemonchiffon", 
	"palegreen", "lightskyblue", "darksalmon", "plum"])
fig.tight_layout()
fig.show()
"""
